// Copyright (c) 2018-2026  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

#ifndef TERRA_TBB_HELPER_H
#define TERRA_TBB_HELPER_H

#if defined(HAVE_TBB) && !defined(USE_TBB)
#define USE_TBB
#endif

#if defined(USE_TBB)
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_arena.h>
#endif

#include "spatBase.h"

// Single switch + thread-cap for TBB-parallel kernels.
//   opt.parallel  : on/off (also gated by the caller; this header does not
//                   re-check, callers wrap calls in `if (opt.parallel)`)
//   opt.threads   : maximum number of threads. 0 = no cap (TBB default).
//
// Use `terra_parallel_for(opt, range, body)` instead of calling
// tbb::parallel_for directly so that opt.threads is honored consistently.
#if defined(USE_TBB)

template <typename Range, typename Body>
inline void terra_parallel_for(const SpatOptions &opt,
                               const Range &range, const Body &body) {
	if (opt.threads > 0) {
		tbb::task_arena arena((int) opt.threads);
		arena.execute([&]{ tbb::parallel_for(range, body); });
	} else {
		tbb::parallel_for(range, body);
	}
}

#endif // USE_TBB

#endif // TERRA_TBB_HELPER_H
