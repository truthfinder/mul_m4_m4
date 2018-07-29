#pragma once
#ifndef timer_h__
#define timer_h__

//std::chrono::duration_cast<std::chrono::seconds>(tp-tp).count()

template <typename T> class check_ts { // temporary class for performance debugging purposes
public:
	static __int64 freq() noexcept { return T::freq(); }

	check_ts(bool const started = true) noexcept : t(0ll), r(0ll), f(1./T::freq()) { if (started) { start(); } }
	~check_ts() noexcept { stop(); }

	__int64      start() noexcept { return t = T::tics_start(); }
	__int64       stop() noexcept { return r = T::tics_stop() - t; }
	double    stop_sec() noexcept { return f * stop(); }
	__int64      reset() noexcept { __int64 c = T::tics_start(); r = c - t; t = c; return r; }
	double   reset_sec() noexcept { return f * reset(); }
	__int64     passed() const noexcept { return r; }
	__int64  passed_ms() const noexcept {
		std::chrono::duration<__int64, std::ratio<1, T::freq()>> const u(passed());
		return std::chrono::duration_cast<std::chrono::milliseconds>(u).count();
	}
	__int64 passed_mks() const noexcept {
		std::chrono::duration<__int64, std::ratio<1, T::freq()>> const u(passed());
		return std::chrono::duration_cast<std::chrono::microseconds>(u).count();
	}
	__int64 passed_ns() const noexcept {
		return passed();
	}
	double  passed_sec() const noexcept { return f * passed(); }
	__int64    current() const noexcept { return T::tics() - t; }
	double current_sec() const noexcept { return f * current(); }

	check_ts& operator += (check_ts const& ts) noexcept { r += ts.r; return *this; }
	check_ts& operator -= (check_ts const& ts) noexcept { r -= ts.r; return *this; }
	check_ts& operator += ( __int64 const   d) noexcept { r +=    d; return *this; }
	check_ts& operator -= ( __int64 const   d) noexcept { r -=    d; return *this; }

private:
	__int64 t, r;
	double f;
};

class hrc_timer { // temporary class for performance debugging purposes
	using clock = std::chrono::high_resolution_clock;
	//using frq = std::chrono::nanoseconds = clock::type?;
public:
	static __int64 tics_start() noexcept { return clock::now().time_since_epoch().count(); }
	static __int64 tics_stop() noexcept { return tics_start(); }
	static constexpr __int64 freq() noexcept { return clock::period::den; }
};

class rdtsc_timer {
public:
	static __int64 tics_start() noexcept {
		int info[4];
		__cpuid(info, 1);
		return __rdtsc();
	}
	static __int64 tics_stop() noexcept {
		int info[4];
		unsigned int aux;
		__int64 r = __rdtscp(&aux);
		__cpuid(info, 1);
		return r;
	}
	static __int64 freq() noexcept {
		return 3'700'000'000ll;
	}
};

typedef check_ts<hrc_timer> hts;
typedef check_ts<rdtsc_timer> rts;

#endif // timer_h__
