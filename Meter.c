/*
htop - Meter.c
(C) 2004-2011 Hisham H. Muhammad
Released under the GNU GPLv2+, see the COPYING file
in the source distribution for its full text.
*/

#include "config.h" // IWYU pragma: keep

#include "Meter.h"

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CRT.h"
#include "Macros.h"
#include "Object.h"
#include "ProvideCurses.h"
#include "RichString.h"
#include "Settings.h"
#include "XUtils.h"


#define GRAPH_HEIGHT 4 /* Unit: rows (lines) */

const MeterClass Meter_class = {
   .super = {
      .extends = Class(Object)
   }
};

Meter* Meter_new(const Machine* host, unsigned int param, const MeterClass* type) {
   Meter* this = xCalloc(1, sizeof(Meter));
   Object_setClass(this, type);
   this->h = 1;
   this->param = param;
   this->host = host;
   this->curItems = type->maxItems;
   this->curAttributes = NULL;
   this->values = type->maxItems ? xCalloc(type->maxItems, sizeof(double)) : NULL;
   this->total = type->total;
   this->caption = xStrdup(type->caption);
   if (Meter_initFn(this)) {
      Meter_init(this);
   }
   Meter_setMode(this, type->defaultMode);
   return this;
}

static const char* const Meter_prefixes = "KMGTPEZY";

int Meter_humanUnit(char* buffer, unsigned long int value, size_t size) {
   const char* prefix = Meter_prefixes;
   unsigned long int powi = 1;
   unsigned int powj = 1, precision = 2;

   for (;;) {
      if (value / 1024 < powi)
         break;

      if (prefix[1] == '\0')
         break;

      powi *= 1024;
      ++prefix;
   }

   if (*prefix == 'K')
      precision = 0;

   for (; precision > 0; precision--) {
      powj *= 10;
      if (value / powi < powj)
         break;
   }

   return snprintf(buffer, size, "%.*f%c", precision, (double) value / powi, *prefix);
}

static void GraphMeterMode_printScale(int scaleExp) {
   if (scaleExp > 86) { // > 99 yotta
      addstr("inf");
   } else if (scaleExp < 10) { // <= 512
      if (scaleExp < 0) {
         // This is unimplemented because GraphMeterMode_generateRecord sets
         // scaleExp with a minimum value of 0.
         return;
      }
      printw("%3u", 1U << scaleExp);
   } else if (scaleExp % 10 <= 6) { // {1|2|4|8|16|32|64}{K|M|G|T|P|E|Z|Y}
      printw("%2u%c", 1U << (scaleExp % 10), Meter_prefixes[(scaleExp / 10) - 1]);
   } else {
      // To keep the 3 character limit, we express in the "fraction" format
      // such as "M/8" (=128K), "G/4" (=256M) and "T/2" (=512G).
      printw("%c/%u", Meter_prefixes[scaleExp / 10], 1U << (10 - (scaleExp % 10)));
   }
}

void Meter_delete(Object* cast) {
   if (!cast)
      return;

   Meter* this = (Meter*) cast;
   if (Meter_doneFn(this)) {
      Meter_done(this);
   }
   if (this->drawData)
      free(this->drawData->cells);
   free(this->drawData);
   free(this->caption);
   free(this->values);
   free(this);
}

void Meter_setCaption(Meter* this, const char* caption) {
   free_and_xStrdup(&this->caption, caption);
}

static inline void Meter_displayBuffer(const Meter* this, RichString* out) {
   if (Object_displayFn(this)) {
      Object_display(this, out);
   } else {
      RichString_writeWide(out, CRT_colors[Meter_attributes(this)[0]], this->txtBuffer);
   }
}

void Meter_setMode(Meter* this, int modeIndex) {
   if (modeIndex > 0 && modeIndex == this->mode) {
      return;
   }

   if (!modeIndex) {
      modeIndex = 1;
   }

   assert(modeIndex < LAST_METERMODE);
   if (Meter_defaultMode(this) == CUSTOM_METERMODE) {
      this->draw = Meter_drawFn(this);
      if (Meter_updateModeFn(this)) {
         Meter_updateMode(this, modeIndex);
      }
   } else {
      assert(modeIndex >= 1);
      free(this->drawData);
      this->drawData = NULL;

      const MeterMode* mode = Meter_modes[modeIndex];
      this->draw = mode->draw;
      this->h = mode->h;
   }
   this->mode = modeIndex;
}

ListItem* Meter_toListItem(const Meter* this, bool moving) {
   char mode[20];
   if (this->mode) {
      xSnprintf(mode, sizeof(mode), " [%s]", Meter_modes[this->mode]->uiName);
   } else {
      mode[0] = '\0';
   }
   char name[32];
   if (Meter_getUiNameFn(this))
      Meter_getUiName(this, name, sizeof(name));
   else
      xSnprintf(name, sizeof(name), "%s", Meter_uiName(this));
   char buffer[50];
   xSnprintf(buffer, sizeof(buffer), "%s%s", name, mode);
   ListItem* li = ListItem_new(buffer, 0);
   li->moving = moving;
   return li;
}

static double Meter_computeSum(const Meter* this) {
   double sum = 0.0;
   for (uint8_t i = 0; i < this->curItems; i++) {
      if (isgreater(this->values[i], 0.0)) {
         sum += this->values[i];
         // Cap the sum to a finite value so we can compute percentages.
         if (!(sum <= DBL_MAX)) {
            return DBL_MAX;
         }
      }
   }
   return sum;
}

/* ---------- TextMeterMode ---------- */

static void TextMeterMode_draw(Meter* this, int x, int y, int w) {
   const char* caption = Meter_getCaption(this);
   attrset(CRT_colors[METER_TEXT]);
   mvaddnstr(y, x, caption, w);
   attrset(CRT_colors[RESET_COLOR]);

   int captionLen = strlen(caption);
   x += captionLen;
   w -= captionLen;
   if (w <= 0)
      return;

   RichString_begin(out);
   Meter_displayBuffer(this, &out);
   RichString_printoffnVal(out, y, x, 0, w);
   RichString_delete(&out);
}

/* ---------- BarMeterMode ---------- */

static const char BarMeterMode_characters[] = "|#*@$%&.";

static void BarMeterMode_draw(Meter* this, int x, int y, int w) {
   const char* caption = Meter_getCaption(this);
   attrset(CRT_colors[METER_TEXT]);
   int captionLen = 3;
   mvaddnstr(y, x, caption, captionLen);
   x += captionLen;
   w -= captionLen;
   attrset(CRT_colors[BAR_BORDER]);
   mvaddch(y, x, '[');
   w--;
   mvaddch(y, x + MAXIMUM(w, 0), ']');
   w--;
   attrset(CRT_colors[RESET_COLOR]);

   x++;

   if (w < 1)
      return;

   // The text in the bar is right aligned;
   // Pad with maximal spaces and then calculate needed starting position offset
   RichString_begin(bar);
   RichString_appendChr(&bar, 0, ' ', w);
   RichString_appendWide(&bar, 0, this->txtBuffer);
   int startPos = RichString_sizeVal(bar) - w;
   if (startPos > w) {
      // Text is too large for bar
      // Truncate meter text at a space character
      for (int pos = 2 * w; pos > w; pos--) {
         if (RichString_getCharVal(bar, pos) == ' ') {
            while (pos > w && RichString_getCharVal(bar, pos - 1) == ' ')
               pos--;
            startPos = pos - w;
            break;
         }
      }

      // If still too large, print the start not the end
      startPos = MINIMUM(startPos, w);
   }
   assert(startPos >= 0);
   assert(startPos <= w);
   assert(startPos + w <= RichString_sizeVal(bar));

   int blockSizes[10];
   int blockSizeSum = 0;

   // First draw in the bar[] buffer...
   int offset = 0;
   double total = this->total;
   if (total <= 0.0) {
      // Update "total" value for dynamic scale.
      // "total" is only updated in Bar meter mode.
      double sum = Meter_computeSum(this);
      if (sum > -this->total)
         this->total = -sum;

      total = -this->total;
   }
   for (uint8_t i = 0; i < this->curItems; i++) {
      double value = this->values[i];
      value = isgreater(value, 0.0) ? MINIMUM(value, total) : 0.0;
      if (value > 0) {
         blockSizes[i] = ceil((value / total) * w);
      } else {
         blockSizes[i] = 0;
      }

      if (Meter_comprisedValues(this)) {
         blockSizes[i] = MAXIMUM(blockSizes[i] - blockSizeSum, 0);
         blockSizeSum += blockSizes[i];
      }

      int nextOffset = offset + blockSizes[i];
      // (Control against invalid values)
      nextOffset = CLAMP(nextOffset, 0, w);
      for (int j = offset; j < nextOffset; j++)
         if (RichString_getCharVal(bar, startPos + j) == ' ') {
            if (CRT_colorScheme == COLORSCHEME_MONOCHROME) {
               assert(i < strlen(BarMeterMode_characters));
               RichString_setChar(&bar, startPos + j, BarMeterMode_characters[i]);
            } else {
               RichString_setChar(&bar, startPos + j, '|');
            }
         }
      offset = nextOffset;
   }

   // ...then print the buffer.
   offset = 0;
   for (uint8_t i = 0; i < this->curItems; i++) {
      int attr = this->curAttributes ? this->curAttributes[i] : Meter_attributes(this)[i];
      RichString_setAttrn(&bar, CRT_colors[attr], startPos + offset, blockSizes[i]);
      RichString_printoffnVal(bar, y, x + offset, startPos + offset, MINIMUM(blockSizes[i], w - offset));
      offset += blockSizes[i];
      offset = CLAMP(offset, 0, w);
   }
   if (offset < w) {
      RichString_setAttrn(&bar, CRT_colors[BAR_SHADOW], startPos + offset, w - offset);
      RichString_printoffnVal(bar, y, x + offset, startPos + offset, w - offset);
   }

   RichString_delete(&bar);

   move(y, x + w + 1);
   attrset(CRT_colors[RESET_COLOR]);
}

/* ---------- GraphMeterMode ---------- */

static void GraphMeterMode_createGraphData(Meter* this, unsigned int numRecords) {
   assert(this->drawData == NULL);

   GraphData* data = xCalloc(1, sizeof(GraphData));
   data->numRecords = numRecords;

   bool isOneColorGraph = GRAPH_HEIGHT >= 2 && Meter_maxItems(this) <= 1;

   // recordNumCells determines number of cells needed to allocate for a graph
   // data record. A record corresponds to a terminal column when drawn.
   size_t recordNumCells;
   if (isOneColorGraph) {
      // Index 0: log2 of scale
      // Index 1: Number of dots total
      recordNumCells = 2;
   } else if (this->total > 0.0) {
      // Percent graph
      // Index 0..GRAPH_HEIGHT: Color cells
      recordNumCells = GRAPH_HEIGHT;
   } else {
      // Non-percentage & dynamic scale
      // Index 0: log2 of scale
      // Index 1..(2*GRAPH_HEIGHT-1): Color cells
      recordNumCells = 2 * GRAPH_HEIGHT;
   }
   assert(recordNumCells <= UINT16_MAX);
   data->recordNumCells = (uint16_t)recordNumCells;

   GraphColorCell* cells = xCalloc(numRecords * recordNumCells, sizeof(GraphColorCell));
   if (!isOneColorGraph) {
      for (size_t i = 0; i < numRecords * recordNumCells; i++) {
         if (this->total > 0.0 || i % recordNumCells > 0)
            cells[i].cell.itemIndex = UINT8_MAX;
      }
   }
   data->cells = cells;
   this->drawData = data;
}

static uint8_t GraphMeterMode_findTopCellItem(Meter* this, double scaledTotal, int topCell) {
   double sum = 0.0;
   double bottom = 0.0;
   double maxArea = 0.0;
   uint8_t topCellItem = this->curItems - 1;
   for (uint8_t i = 0; i < this->curItems && sum < DBL_MAX; i++) {
      if (!isgreater(this->values[i], 0.0))
         continue;

      sum += this->values[i];
      if (!(sum <= DBL_MAX))
         sum = DBL_MAX;

      double top = (sum / scaledTotal) * GRAPH_HEIGHT;
      if (top > topCell) {
         // Find the item that occupies the largest area of the top cell.
         // Favor item with higher index in case of a tie.
         if (top - bottom >= maxArea) {
            topCellItem = i;
            maxArea = top > 0.0 ? top - bottom : top - topCell;
         }
      }
      bottom = top;
   }
   return topCellItem;
}

static void GraphMeterMode_computeColors(Meter* this, GraphColorCell* record, int numDots, double scaledTotal, int scaleExp, bool isPercentGraph) {
   assert(numDots > 0 && numDots <= (GRAPH_HEIGHT * 8));

   // If there is a "top cell" which won't be completely filled, determine
   // its color first.
   int topCell = (numDots - 1) / 8;
   uint8_t dotAlignment = 2;
   unsigned int blanksAtTopCell = (8 + topCell * 8 - numDots) / dotAlignment * dotAlignment;

   bool hasPartialTopCell = false;
   if (blanksAtTopCell > 0) {
      hasPartialTopCell = true;
   } else if (!isPercentGraph && topCell == ((GRAPH_HEIGHT - 1) >> scaleExp) && topCell % 2 == 0) {
      // This "top cell" is rendered as full in one scale, but partial in the
      // next scale. (Only happens when GRAPH_HEIGHT is not a power of two.)
      hasPartialTopCell = true;
   }

   double topCellArea = 0.0;
   uint8_t topCellItem = this->curItems - 1;
   if (hasPartialTopCell) {
      topCellArea = (8 - blanksAtTopCell) / 8.0;
      topCellItem = GraphMeterMode_findTopCellItem(this, scaledTotal, topCell);
   }
   topCell += 1; // This index points to a cell that would be blank.

   // Compute colors of the rest of the cells, using the largest remainder
   // method (a.k.a. Hamilton's method).
   // The Hare quota is (scaledTotal / GRAPH_HEIGHT).
   int paintedHigh = topCell + (int)topCellItem + 1;
   int paintedLow = 0;
   double threshold = 0.5;
   double thresholdHigh = 1.0;
   double thresholdLow = 0.0;
   // Tiebreak 1: Favor items with less number of cells. (Top cell is not
   // included in the count.)
   int cellLimit = topCell;
   int cellLimitHigh = topCell;
   int cellLimitLow = 0.0;
   // Tiebreak 2: Favor items whose indices are lower.
   uint8_t tiedItemLimit = topCellItem;
   while (true) {
      double sum = 0.0;
      double bottom = 0.0;
      int cellsPainted = 0;
      double nextThresholdHigh = 0.0;
      double nextThresholdLow = 1.0;
      int nextCellLimitHigh = 0;
      int nextCellLimitLow = topCell;
      uint8_t numTiedItems = 0;
      for (uint8_t i = 0; i <= topCellItem && sum < DBL_MAX; i++) {
         if (!isgreater(this->values[i], 0.0))
            continue;

         sum += this->values[i];
         if (!(sum <= DBL_MAX))
            sum = DBL_MAX;

         double top = (sum / scaledTotal) * GRAPH_HEIGHT;
         double area = top - bottom;
         double remainder = area;
         if (i == topCellItem) {
            remainder -= topCellArea;
            if (!(remainder >= 0.0))
               remainder = 0.0;
         }
         int numCells = (int)remainder;
         remainder -= numCells;

         // Whether the item will receive an extra cell or be truncated
         if (remainder >= threshold) {
            if (remainder > threshold) {
               if (remainder < nextThresholdLow)
                  nextThresholdLow = remainder;
               numCells++;
            } else if (numCells <= cellLimit) {
               if (numCells < cellLimit) {
                  if (numCells > nextCellLimitHigh)
                     nextCellLimitHigh = numCells;
                  numCells++;
               } else {
                  numTiedItems++;
                  if (numTiedItems <= tiedItemLimit) {
                     numCells++;
                  } else {
                     remainder = 0.0;
                  }
               }
            } else {
               if (numCells < nextCellLimitLow)
                  nextCellLimitLow = numCells;
               remainder = 0.0;
            }
         } else {
            if (remainder > nextThresholdHigh)
               nextThresholdHigh = remainder;
            remainder = 0.0;
         }

         // Paint cells to the buffer
         uint8_t blanksAtEnd = 0;
         if (i == topCellItem && topCellArea > 0.0) {
            numCells++;
            if (area < topCellArea) {
               remainder = MAXIMUM(area, 0.25);
               blanksAtEnd = (uint8_t)blanksAtTopCell;
            }
         } else if (cellsPainted + numCells >= topCell) {
            blanksAtEnd = 0;
         } else if (cellsPainted <= 0 || bottom <= cellsPainted) {
            blanksAtEnd = ((uint8_t)((1.0 - remainder) * 8.0) % 8);
         } else if (cellsPainted + numCells > top) {
            assert(cellsPainted + numCells - top < 1.0);
            blanksAtEnd = (uint8_t)((cellsPainted + numCells - top) * 8.0);
         }

         unsigned int blanksAtStart = 0;
         if (cellsPainted > 0) {
            blanksAtStart = ((uint8_t)((1.0 - remainder) * 8.0) % 8 - blanksAtEnd);
         }

         while (numCells > 0 && cellsPainted < topCell) {
            unsigned int offset = (unsigned int)cellsPainted;
            if (!isPercentGraph) {
               offset = (offset * 2 + 1) << scaleExp;
            }

            record[offset].cell.itemIndex = (uint8_t)i;
            record[offset].cell.details = 0xFF;

            if (blanksAtStart > 0) {
               assert(blanksAtStart < 8);
               record[offset].cell.details >>= blanksAtStart;
               blanksAtStart = 0;
            }

            if (cellsPainted == topCell - 1) {
               assert(blanksAtTopCell < 8);
               record[offset].cell.details &= 0xFF << blanksAtTopCell;
            } else if (numCells == 1) {
               assert(blanksAtEnd < 8);
               record[offset].cell.details &= 0xFF << blanksAtEnd;
            }

            numCells--;
            cellsPainted++;
         }
         cellsPainted += numCells;

         bottom = top;
      }

      if (cellsPainted == topCell)
         break;

      // Set new bounds and threshold
      if (cellsPainted > topCell) {
         paintedHigh = cellsPainted;
         if (thresholdLow >= thresholdHigh) {
            if (cellLimitLow >= cellLimitHigh) {
               assert(tiedItemLimit >= topCellItem);
               tiedItemLimit = numTiedItems - (uint8_t)(cellsPainted - topCell);
            } else {
               assert(cellLimitHigh > cellLimit);
               cellLimitHigh = cellLimit;
            }
         } else {
            assert(thresholdLow < threshold);
            thresholdLow = threshold;
         }
      } else {
         paintedLow = cellsPainted + 1;
         if (thresholdLow >= thresholdHigh) {
            assert(cellLimitLow < cellLimitHigh);
            assert(cellLimitLow < nextCellLimitLow);
            cellLimitLow = nextCellLimitLow;
            nextCellLimitHigh = cellLimitHigh;
         } else {
            assert(cellLimit >= topCell);
            assert(thresholdHigh > nextThresholdHigh);
            thresholdHigh = nextThresholdHigh;
            nextThresholdLow = thresholdLow;
         }
      }
      threshold = thresholdHigh;
      if (thresholdLow >= thresholdHigh) {
         cellLimit = cellLimitLow;
         if (cellLimitLow < cellLimitHigh && paintedHigh > paintedLow) {
            cellLimit += ((cellLimitHigh - cellLimitLow) *
                         (topCell - paintedLow) /
                         (paintedHigh - paintedLow));
            if (cellLimit > nextCellLimitHigh)
               cellLimit = nextCellLimitHigh;
         }
      } else {
         if (paintedHigh > paintedLow) {
            threshold -= ((thresholdHigh - thresholdLow) *
                         (topCell - paintedLow) /
                         (paintedHigh - paintedLow));
            if (threshold < nextThresholdLow)
               threshold = nextThresholdLow;
         }
      }
      assert(threshold <= thresholdHigh);
   }
}

static void GraphMeterMode_generateRecord(Meter* this) {
   assert(this->drawData != NULL);

   // Move previous records
   GraphData* data = this->drawData;
   size_t numRecords = data->numRecords;
   size_t recordNumCells = data->recordNumCells;
   memmove(&data->cells[0],
           &data->cells[1 * recordNumCells],
           ((numRecords - 1) * recordNumCells * sizeof(GraphColorCell)));
   GraphColorCell* record = &data->cells[(numRecords - 1) * recordNumCells];

   // We use braille dots to depict values of each data item. We can only draw
   // one color per character (i.e. cell).

   // Compute "sum" and "total"
   bool isPercentGraph = this->total > 0.0;
   double sum = Meter_computeSum(this);
   assert(sum >= 0.0 && sum <= DBL_MAX);
   double total;
   int exp = 0;
   if (isPercentGraph) {
      total = MAXIMUM(this->total, sum);
   } else {
      (void) frexp(sum, &exp);
      if (exp < 0) {
         exp = 0;
      } else {
         // Overflow case (non-IEEE float; mentioned for completeness)
         // 2^32768 is approximately 1.41546 * 10^9864.
         if (DBL_MAX_10_EXP >= 9864 && exp > INT16_MAX) {
            exp = INT16_MAX;
         }
      }
      assert(exp <= INT16_MAX);
      record[0].scaleExp = (int16_t)exp;
      total = ldexp(1.0, exp);
   }
   if (!(total <= DBL_MAX))
      total = DBL_MAX; // Either this->total or ldexp(1.0, exp) overflows.

   int numDots = (int) ceil((sum / total) * (GRAPH_HEIGHT * 8));
   if (sum > 0.0 && numDots <= 0)
      numDots = 1; // (sum / total) underflow

   bool isOneColorGraph = GRAPH_HEIGHT >= 2 && Meter_maxItems(this) <= 1;
   if (isOneColorGraph) {
      assert(GRAPH_HEIGHT <= UINT16_MAX / 8);
      assert(recordNumCells == 2);
      assert(record[0].scaleExp == exp);
      assert(numDots <= UINT16_MAX);
      record[1].numDots = (uint16_t)numDots;
      return;
   }

   // Clear cells
   size_t offset = (size_t)((numDots - 1) / 8 + 1); // Round up
   assert(offset <= GRAPH_HEIGHT);
   if (!isPercentGraph) {
      assert(recordNumCells == 2 * GRAPH_HEIGHT);
      offset = offset * 2 + 1;
   }
   for (; offset < recordNumCells; offset++) {
      record[offset].cell.itemIndex = UINT8_MAX;
      record[offset].cell.details = 0x00;
   }
   if (sum <= 0.0)
      return;

   int scaleExp = 0;
   double scaledTotal;
   do {
      scaledTotal = ldexp(total, scaleExp);
      if (!(scaledTotal <= DBL_MAX))
         scaledTotal = DBL_MAX;
      numDots = (int) ceil((sum / scaledTotal) * (GRAPH_HEIGHT * 8));
      if (numDots <= 0)
         numDots = 1; // (sum / scaledTotal) underflow

      GraphMeterMode_computeColors(this, record, numDots, scaledTotal, scaleExp, isPercentGraph);
      scaleExp++;
   } while (!isPercentGraph && scaledTotal < DBL_MAX && (1U << scaleExp) < GRAPH_HEIGHT * 2);
}

static uint8_t GraphMeterMode_scaleCellDetails(uint8_t details, int scaleExp) {
   if (scaleExp < 1) {
      assert(scaleExp == 0);
      // No scaling
      return details;
   }
   if (scaleExp == 1 && (details & 0x0F) != 0x00) {
      // Shrink the cell to half height (bits 0 to 3 are zero).
      // For readability, bits 4 and 5 are always set.
      uint8_t newDetails = 0x30;
      if (popCount8(details) >= 5)
         newDetails |= 0x40; // Bit 6
      if (details >= 0x7F)
         newDetails |= 0x80; // Bit 7
      return newDetails;
   }
   // Data will occupy a quarter of the cell or less (bits 0 to 5 are zero).
   // If the cell is non-zero, show at least a quarter (bits 6 and 7 set).
   if (details != 0x00)
      return 0xC0;
   return 0x00;
}

static int GraphMeterMode_lookupCell(const Meter* this, unsigned int recordIndex, unsigned int h, int scaleExp, uint8_t* details) {
   assert(this != NULL);
   assert(details != NULL);

   // Empty cell
   *details = (h == 0) ? 0xC0 : 0x00;

   uint8_t maxItems = Meter_maxItems(this);
   if (maxItems <= 0)
      return BAR_SHADOW;

   assert(this->drawData != NULL);
   GraphData* data = this->drawData;
   assert(recordIndex < data->numRecords);
   uint16_t recordNumCells = data->recordNumCells;
   const GraphColorCell* record = &data->cells[recordIndex * recordNumCells];
   const GraphColorCell* colorCell;
   uint8_t itemIndex;
   if (GRAPH_HEIGHT >= 2 && maxItems <= 1) {
      // One-color graph
      assert(recordNumCells == 2);

      uint16_t numDots = record[1].numDots;
      if (numDots <= 0)
         return BAR_SHADOW;

      int expDifference = scaleExp - record[0].scaleExp;
      assert(expDifference >= 0);

      // Scale according to exponent difference. Round up.
      if (expDifference >= 16) {
         numDots = 1;
      } else {
         numDots = ((numDots - 1) >> expDifference) + 1;
      }
      if (h * 8 >= numDots)
         return BAR_SHADOW;

      *details = 0xFF;
      if ((h + 1) * 8 > numDots) {
         uint8_t dotAlignment = 2;
         unsigned int blanksAtTopCell = ((h + 1) * 8 - numDots) / dotAlignment * dotAlignment;
         *details <<= blanksAtTopCell;
      }
      itemIndex = 0;
   } else if (this->total > 0.0) {
      // Percent graph
      assert(recordNumCells == GRAPH_HEIGHT);

      colorCell = &record[h];
      *details = colorCell->cell.details;
      itemIndex = colorCell->cell.itemIndex;
   } else {
      // Dynamic scale: A record may be rendered in different scales depending
      // on the largest value of a given time frame. The colors need to be
      // determined in different sizes for the same record.
      // An array indexing trick is used in storing cells' color data for
      // different scales together.
      //
      // Example: GRAPH_HEIGHT = 6, n = sum value rounded up to a power of 2.
      //
      //   scale   1*n  2*n  4*n  8*n 16*n
      // ---------------------------------
      //   array  [11]
      // indices   [9]
      //           [7]                      (Cells with array indices greater
      //           [5] [10]                 than (2*GRAPH_HEIGHT-1) are
      //           [3]  [6] (12)            computed from cells of a lower
      //           [1]  [2]  [4]  [8] (16)  scale. They are not stored.)
      assert(recordNumCells == 2 * GRAPH_HEIGHT);

      int expDifference = scaleExp - record[0].scaleExp;
      assert(expDifference >= 0);
      if (h > ((unsigned int)(GRAPH_HEIGHT - 1) >> expDifference))
         return BAR_SHADOW;
      unsigned int a = (h << (expDifference + 1));
      if (h == ((unsigned int)(GRAPH_HEIGHT - 1) >> expDifference)) {
         colorCell = &record[a + powerOf2Floor(2 * GRAPH_HEIGHT - 1 - a)];
         *details = GraphMeterMode_scaleCellDetails(colorCell->cell.details, expDifference);
      } else {
         colorCell = &record[a + (1U << expDifference)];
         *details = colorCell->cell.details;
      }
      itemIndex = colorCell->cell.itemIndex;
   }

   if (h == 0)
      *details |= 0xC0;
   if (itemIndex == UINT8_MAX)
      return BAR_SHADOW;
   assert(itemIndex < maxItems);
   return Meter_attributes(this)[itemIndex];
}

static char GraphMeterMode_detailsToDotsAscii(uint8_t details) {
   // Use ':' if the total number of dots is 5 or more.
   if (popCount8(details) >= 5)
      return ':';
   // Determine which half has more dots than the other.
   uint8_t inverted = details ^ 0x0F;
   int difference = popCount8(inverted) - 4;
   if (difference < 0)
      return '\'';
   if (difference > 0)
      return '.';
   // Upper and lower half have equal number of dots. Give favor to dots more
   // distanced away from the center of the cell.
   // Reverse bits 0 to 3 and subtract it from bits 4 to 7. The algorithm is
   // optimized for IA-32.
   uint32_t n = (uint32_t)(inverted * 0x010101UL);
   n = (uint32_t)((n & 0xF20508UL) * 0x028821UL);
   difference = ((n >> 20) & 0x1F) - 0x0F;
   if (difference < 0)
      return '\'';
   return '.';
}

static void GraphMeterMode_printDots(uint8_t details) {
   if (details == 0x00) {
      // Use ASCII space instead. A braille blank character may display as a
      // substitute block ("tofu").
      addch(' ');
      return;
   }
#ifdef HAVE_LIBNCURSESW
   if (CRT_utf8) {
      // Convert GraphColorCell.details bit representation to Unicode braille
      // dot ordering.
      //   (Bit0) a b (Bit3)  From:        h g f e d c b a (binary)
      //   (Bit1) c d (Bit4)               | | |  X   X  |
      //   (Bit2) e f (Bit5)               | | | | \ / | |
      //   (Bit6) g h (Bit7)               | | | |  X  | |
      //                      To: 0x2800 + h g f d b e c a
      //
      // Braille patterns (U+2800 - U+28FF) in UTF-8: (E2 A0 80 - E2 A3 BF)
      char sequence[4] = "\xE2\xA0\x80";
      // Bits 6 and 7 is in the second byte of the UTF-8 sequence.
      sequence[1] |= details >> 6;
      // Bits 0 to 5 in the third byte. The algorithm is optimized for IA-32
      // and is smaller than naive bit extraction and masking code.
      //
      // Make copies of same data.
      uint32_t n = (uint32_t)(details * 0x01010101UL);
      // Mask bits. (The mask value and multiplier must be chosen so that there
      // are no bit overlap during multiplication.)
      n &= 0x08211204UL;
      // Shift and add.
      n *= 0x00422081UL;
      sequence[2] |= (n >> 23) & 0x3F;
      addstr(sequence);
      return;
   }
#endif
   // ASCII mode
   addch(GraphMeterMode_detailsToDotsAscii(details));
}

static void GraphMeterMode_draw(Meter* this, int x, int y, int w) {
   const int captionLen = 3;
   const char* caption = Meter_getCaption(this);
   attrset(CRT_colors[METER_TEXT]);
   mvaddnstr(y, x, caption, captionLen);

   uint8_t maxItems = Meter_maxItems(this);
   if (maxItems <= 0)
      return;

   if (this->drawData == NULL) {
      GraphMeterMode_createGraphData(this, METER_GRAPHDATA_SIZE);
   }
   GraphData* data = this->drawData;

   const Machine* host = this->host;
   if (!timercmp(&host->realtime, &(data->time), <)) {
      int globalDelay = host->settings->delay;
      struct timeval delay = { .tv_sec = globalDelay / 10, .tv_usec = (globalDelay % 10) * 100000L };
      timeradd(&host->realtime, &delay, &(data->time));

      GraphMeterMode_generateRecord(this);
   }

   // How many values (and columns) we can draw for this graph
   w -= captionLen;
   if (w <= 0)
      return;
   unsigned int index = 0;
   int col = 0;
   if (data->numRecords >= (unsigned int)w) {
      index = data->numRecords - (unsigned int)w;
   } else {
      col = w - (int)data->numRecords;
   }

   // If it's not percent graph, determine the scale
   int exp = 0;
   bool isPercentGraph = this->total > 0.0;
   if (!isPercentGraph) {
      uint16_t recordNumCells = data->recordNumCells;
      for (unsigned int j = index; j < data->numRecords; j++) {
         const GraphColorCell* record = &data->cells[j * recordNumCells];
         if (record[0].scaleExp > exp) {
            exp = record[0].scaleExp;
         }
      }
   }

   // Print scale
   if (GRAPH_HEIGHT >= 3) {
      move(y + 1, x);
      if (isPercentGraph) {
         addstr("  %");
      } else {
         GraphMeterMode_printScale(exp);
      }
   }
   x += captionLen;

   // Print the graph
   x += col;
   for (int line = 0; line < GRAPH_HEIGHT; line++) {
      for (col = 0; index + (unsigned int)col < data->numRecords; col++) {
         uint8_t details;
         int colorIdx = GraphMeterMode_lookupCell(this, index + (unsigned int)col, GRAPH_HEIGHT - 1 - line, exp, &details);
         move(y + line, x + col);
         attrset(CRT_colors[colorIdx]);
         GraphMeterMode_printDots(details);
      }
   }
   attrset(CRT_colors[RESET_COLOR]);
}

/* ---------- LEDMeterMode ---------- */

static const char* const LEDMeterMode_digitsAscii[] = {
   " __ ", "    ", " __ ", " __ ", "    ", " __ ", " __ ", " __ ", " __ ", " __ ",
   "|  |", "   |", " __|", " __|", "|__|", "|__ ", "|__ ", "   |", "|__|", "|__|",
   "|__|", "   |", "|__ ", " __|", "   |", " __|", "|__|", "   |", "|__|", " __|"
};

#ifdef HAVE_LIBNCURSESW

static const char* const LEDMeterMode_digitsUtf8[] = {
   "┌──┐", "  ┐ ", "╶──┐", "╶──┐", "╷  ╷", "┌──╴", "┌──╴", "╶──┐", "┌──┐", "┌──┐",
   "│  │", "  │ ", "┌──┘", " ──┤", "└──┤", "└──┐", "├──┐", "   │", "├──┤", "└──┤",
   "└──┘", "  ╵ ", "└──╴", "╶──┘", "   ╵", "╶──┘", "└──┘", "   ╵", "└──┘", " ──┘"
};

#endif

static const char* const* LEDMeterMode_digits;

static void LEDMeterMode_drawDigit(int x, int y, int n) {
   for (int i = 0; i < 3; i++)
      mvaddstr(y + i, x, LEDMeterMode_digits[i * 10 + n]);
}

static void LEDMeterMode_draw(Meter* this, int x, int y, int w) {
#ifdef HAVE_LIBNCURSESW
   if (CRT_utf8)
      LEDMeterMode_digits = LEDMeterMode_digitsUtf8;
   else
#endif
      LEDMeterMode_digits = LEDMeterMode_digitsAscii;

   RichString_begin(out);
   Meter_displayBuffer(this, &out);

   int yText =
#ifdef HAVE_LIBNCURSESW
      CRT_utf8 ? y + 1 :
#endif
      y + 2;
   attrset(CRT_colors[LED_COLOR]);
   const char* caption = Meter_getCaption(this);
   mvaddstr(yText, x, caption);
   int xx = x + strlen(caption);
   int len = RichString_sizeVal(out);
   for (int i = 0; i < len; i++) {
      int c = RichString_getCharVal(out, i);
      if (c >= '0' && c <= '9') {
         if (xx - x + 4 > w)
            break;

         LEDMeterMode_drawDigit(xx, y, c - '0');
         xx += 4;
      } else {
         if (xx - x + 1 > w)
            break;
#ifdef HAVE_LIBNCURSESW
         const cchar_t wc = { .chars = { c, '\0' }, .attr = 0 }; /* use LED_COLOR from attrset() */
         mvadd_wch(yText, xx, &wc);
#else
         mvaddch(yText, xx, c);
#endif
         xx += 1;
      }
   }
   attrset(CRT_colors[RESET_COLOR]);
   RichString_delete(&out);
}

static MeterMode BarMeterMode = {
   .uiName = "Bar",
   .h = 1,
   .draw = BarMeterMode_draw,
};

static MeterMode TextMeterMode = {
   .uiName = "Text",
   .h = 1,
   .draw = TextMeterMode_draw,
};

static MeterMode GraphMeterMode = {
   .uiName = "Graph",
   .h = GRAPH_HEIGHT,
   .draw = GraphMeterMode_draw,
};

static MeterMode LEDMeterMode = {
   .uiName = "LED",
   .h = 3,
   .draw = LEDMeterMode_draw,
};

const MeterMode* const Meter_modes[] = {
   NULL,
   &BarMeterMode,
   &TextMeterMode,
   &GraphMeterMode,
   &LEDMeterMode,
   NULL
};

/* Blank meter */

static void BlankMeter_updateValues(Meter* this) {
   this->txtBuffer[0] = '\0';
}

static void BlankMeter_display(ATTR_UNUSED const Object* cast, ATTR_UNUSED RichString* out) {
}

static const int BlankMeter_attributes[] = {
   DEFAULT_COLOR
};

const MeterClass BlankMeter_class = {
   .super = {
      .extends = Class(Meter),
      .delete = Meter_delete,
      .display = BlankMeter_display,
   },
   .updateValues = BlankMeter_updateValues,
   .defaultMode = TEXT_METERMODE,
   .maxItems = 0,
   .total = 1.0,
   .attributes = BlankMeter_attributes,
   .name = "Blank",
   .uiName = "Blank",
   .caption = ""
};
