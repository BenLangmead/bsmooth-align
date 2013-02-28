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

#ifdef USE_PTHREADS

#include <pthread.h>
#include <iostream>

using namespace std;

void createThread(pthread_t& thr, void *(*func)(void*), void *vp) {
	pthread_attr_t pt_attr;
	pthread_attr_init(&pt_attr);
	pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_JOINABLE);
	pthread_attr_setstacksize(&pt_attr, 2 << 20);
	int ret;
	if((ret = pthread_create(&thr, &pt_attr, func, vp)) != 0) {
		cerr << "Could not create a thread; return from pthread_create: " << ret << endl;
		throw 1;
	}
}

void joinThread(const pthread_t& thr) {
	pthread_join(thr, NULL);
}

#endif
