#pragma once
#ifndef timer_h__
#define timer_h__

//std::chrono::duration_cast<std::chrono::seconds>(tp-tp).count()

template <typename T> class check_ts { // temporary class for performance debugging purposes
public:
	static int64_t freq() noexcept { return T::freq(); }

	check_ts(bool const started = true) noexcept : t(0ll), r(0ll), f(1./T::freq()) { if (started) { start(); } }
	~check_ts() noexcept { stop(); }

	int64_t      start() noexcept { return t = T::tics_start(); }
	int64_t       stop() noexcept { return r = T::tics_stop() - t; }
	double    stop_sec() noexcept { return f * stop(); }
	int64_t      reset() noexcept { int64_t c = T::tics_start(); r = c - t; t = c; return r; }
	double   reset_sec() noexcept { return f * reset(); }
	int64_t     passed() const noexcept { return r; }
	int64_t  passed_ms() const noexcept {
		std::chrono::duration<int64_t, std::ratio<1, T::freq()>> const u(passed());
		return std::chrono::duration_cast<std::chrono::milliseconds>(u).count();
	}
	int64_t passed_mks() const noexcept {
		std::chrono::duration<int64_t, std::ratio<1, T::freq()>> const u(passed());
		return std::chrono::duration_cast<std::chrono::microseconds>(u).count();
	}
	int64_t passed_ns() const noexcept {
		return passed();
	}
	double  passed_sec() const noexcept { return f * passed(); }
	int64_t    current() const noexcept { return T::tics() - t; }
	double current_sec() const noexcept { return f * current(); }

	check_ts& operator += (check_ts const& ts) noexcept { r += ts.r; return *this; }
	check_ts& operator -= (check_ts const& ts) noexcept { r -= ts.r; return *this; }
	check_ts& operator += ( int64_t const   d) noexcept { r +=    d; return *this; }
	check_ts& operator -= ( int64_t const   d) noexcept { r -=    d; return *this; }

private:
	int64_t t, r;
	double f;
};

class hrc_timer { // temporary class for performance debugging purposes
	using clock = std::chrono::high_resolution_clock;
	//using frq = std::chrono::nanoseconds = clock::type?;
public:
	static int64_t tics_start() noexcept { return clock::now().time_since_epoch().count(); }
	static int64_t tics_stop() noexcept { return tics_start(); }
	static constexpr int64_t freq() noexcept { return clock::period::den; }
};

class rdtsc_timer {
public:
	static int64_t tics_start() noexcept {
		int info[4];
#		ifdef WIN32
		__cpuid(info, 1);
#		else
		__cpuid(1, info[0], info[1], info[2], info[3]);
#		endif
		return __rdtsc();
	}
	static int64_t tics_stop() noexcept {
		int info[4];
		unsigned int aux;
		int64_t r = __rdtscp(&aux);
#		ifdef WIN32
		__cpuid(info, 1);
#		else
		__cpuid(1, info[0], info[1], info[2], info[3]);
#		endif
		return r;
	}
	static int64_t freq() noexcept {
		return 3'700'000'000ll;
	}
};

typedef check_ts<hrc_timer> hts;
typedef check_ts<rdtsc_timer> rts;

#endif // timer_h__
