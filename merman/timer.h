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

#ifndef TIMER_H_
#define TIMER_H_

#include <ctime>
#include <iostream>
#include <iomanip>

/**
 * Use time() call to keep track of elapsed time between creation and
 * destruction.  If verbose is true, Timer will print a message showing
 * elapsed time to the given output stream upon destruction.
 */
class Timer {
public:
	Timer(std::ostream& out, const char *msg, bool verbos) :
		_t(time(0)), _out(out), _msg(msg), _verbose(verbos) { }

	/// Optionally print message
	~Timer() {
		if(_verbose) write(_out);
	}

	/// Return elapsed time since Timer object was created
	time_t elapsed() const {
		return time(0) - _t;
	}

	void write(std::ostream& out) {
		time_t passed = elapsed();
		// Print the message supplied at construction time followed
		// by time elapsed formatted HH:MM:SS
		unsigned int hours   = (int)((passed / 60) / 60);
		unsigned int minutes = (int)((passed / 60) % 60);
		unsigned int seconds = (int)(passed % 60);
		out << _msg
			<< std::setfill ('0') << std::setw (2) << hours << ":"
		    << std::setfill ('0') << std::setw (2) << minutes << ":"
		    << std::setfill ('0') << std::setw (2) << seconds << std::endl;
	}

private:
	time_t        _t;
	std::ostream& _out;
	const char   *_msg;
	bool          _verbose;
};

#endif /*TIMER_H_*/
