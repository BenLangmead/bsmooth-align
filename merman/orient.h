/*
 * Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
 *
 * This file is part of Merman.
 *
 * Merman is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Merman is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Merman.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 *  orient.h
 *
 * Consider these scenarios for which "strand" the read may have aligned to:
 *
 * 1. Non-bisulfite case without --genrevcomps:
 *    ========================================
 *
 *  1a. W Case:
 *      ======
 *
 *          offset=0       offset=L
 *          v                     v
 *   W Ref: ----------------------
 *             ||||
 * FW Read:    ----
 *             ^   ^
 *      offset=H   offset=H+D
 *
 *           Reported offset = H
 *  Printed read orientation = FW
 *                 5' end at = H
 *       Printed orientation = +
 * Printed orientation (ext) = W
 *
 *  1b. WR Case:
 *      =======
 *
 *          offset=0       offset=L
 *          v                     v
 *   W Ref: ----------------------
 *             ||||
 * RC Read:    ----
 *             ^   ^
 *      offset=H   offset=H+D
 *
 *           Reported offset = H
 *  Printed read orientation = RC
 *                 5' end at = H + D - 1
 *       Printed orientation = -
 * Printed orientation (ext) = C
 *
 * 2. Non-bisulfite case with --genrevcomps:
 *    =====================================
 *
 *  2a. W Case: Like 1a
 *      ======
 *
 *  2b. C Case:
 *      ======
 *
 *          offset=0       offset=L
 *          v                     v
 *   C Ref: ----------------------
 *             ||||
 * FW Read:    ----
 *             ^   ^
 *      offset=H   offset=H+D
 *
 *           Reported offset = L - (H + D)
 *  Printed read orientation = RC
 *                 5' end at = L - H - 1
 *       Printed orientation = -
 * Printed orientation (ext) = C
 *
 * 3. 2-strand bisulfite case:
 *    =======================
 *
 *  3a. BW Case:
 *      =======
 *
 *          offset=0       offset=L
 *          v                     v
 *  BW Ref: ----------------------
 *             ||||
 * FW Read:    ----
 *             ^   ^
 *      offset=H   offset=H+D
 *
 *           Reported offset = H
 *  Printed read orientation = FW
 *                 5' end at = H
 *       Printed orientation = +
 * Printed orientation (ext) = BW
 *
 *  3b. BC Case:
 *      =======
 *
 *          offset=0       offset=L
 *          v                     v
 *  BC Ref: ----------------------
 *             ||||
 * FW Read:    ----
 *             ^   ^
 *      offset=H   offset=H+D
 *
 *           Reported offset = L - (H + D)
 *  Printed read orientation = RC
 *                 5' end at = L - H - 1
 *       Printed orientation = -
 * Printed orientation (ext) = BC
 *
 * 4. 4-strand bisulfite case:
 *    =======================
 *
 *  4a. BW Case: Like 3a
 *      =======
 *
 *  4b. BC Case: Like 3b
 *      =======
 *
 *  4c. BWR Case:
 *      ========
 *
 *          offset=0       offset=L
 *          v                     v
 *  BW Ref: ----------------------
 *             ||||
 * RC Read:    ----
 *             ^   ^
 *      offset=H   offset=H+D
 *
 *           Reported offset = H
 *  Printed read orientation = RC
 *                 5' end at = H + D - 1
 *       Printed orientation = - (but this is ambiguous)
 * Printed orientation (ext) = BWR
 *
 *  4d. BCR Case:
 *      ========
 *
 *          offset=0       offset=L
 *          v                     v
 *  BC Ref: ----------------------
 *             ||||
 * RC Read:    ----
 *             ^   ^
 *      offset=H   offset=H+D
 *
 *          Reported offset = L - (H + D)
 * Printed read orientation = FW
 *                5' end at = L - (H + D)
 *
 *           Reported offset = L - (H + D)
 *  Printed read orientation = FW
 *                 5' end at = L - (H + D)
 *       Printed orientation = + (but this is ambiguous)
 * Printed orientation (ext) = BCR
 */
 
#include "assert_helpers.h"

#define ORIENT_WATSON    0
#define ORIENT_CRICK     1
#define ORIENT_FW        0
#define ORIENT_RC        2
#define ORIENT_BISULFITE 4

#define ORIENT_WATSON_FW (ORIENT_WATSON | ORIENT_FW)
#define ORIENT_CRICK_FW  (ORIENT_CRICK  | ORIENT_FW)
#define ORIENT_WATSON_RC (ORIENT_WATSON | ORIENT_RC)
#define ORIENT_CRICK_RC  (ORIENT_CRICK  | ORIENT_RC)

extern const char * orientStrs[8];

static inline bool orientWatson(int orient) {
	assert_range(0, 7, orient);
	return ((orient & ORIENT_CRICK) == 0);
}

static inline bool orientCrick(int orient) {
	assert_range(0, 7, orient);
	return !orientWatson(orient);
}

static inline bool orientFw(int orient) {
	assert_range(0, 7, orient);
	return ((orient & ORIENT_RC) == 0);
}

static inline const char *orientStr(int orient) {
	assert_range(0, 7, orient);
	return orientStrs[orient];
}

static inline bool orientRc(int orient) {
	assert_range(0, 7, orient);
	return !orientFw(orient);
}

static inline bool isWatsonFw(int orient) {
	assert_range(0, 7, orient);
	return orient == ORIENT_WATSON_FW;
}

static inline bool isWatsonRc(int orient) {
	assert_range(0, 7, orient);
	return orient == ORIENT_WATSON_RC;
}

static inline bool isCrickFw(int orient) {
	assert_range(0, 7, orient);
	return orient == ORIENT_CRICK_FW;
}

static inline bool isCrickRc(int orient) {
	assert_range(0, 7, orient);
	return orient == ORIENT_CRICK_RC;
}

/**
 * Return true iff the alignment, when lined up against the Watson
 * strand, has its 5p end on the left.
 */
static inline bool is5pLeft(int orient) {
	assert_range(0, 7, orient);
	return (orient == ORIENT_WATSON_FW || orient == ORIENT_CRICK_RC);
}

/**
 * Return false iff an alignment with the given orientation should be printed
 * as the reverse complement of the original read rather than as the original
 * read itself.  Corresponds to cases 1b, 2b, 3b, 4b and 4c in the comment
 * above.
 */
static inline bool orientPrintFw(int orient) {
	assert_range(0, 7, orient);
	return !(isWatsonRc(orient) || isCrickFw(orient));
}
