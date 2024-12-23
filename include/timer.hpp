#pragma once

#include <chrono>
#include <string>
#include <iostream>

namespace AR {
	class Timestep {
	public:
		explicit Timestep(float time = 0.0f) : m_Time(time) {}

		float getMilliseconds() const { return m_Time * 1000.0f; }
		float getSeconds() const { return m_Time; }
		float getFPS() const { return (m_Time > 0.0f) ? (1000.0f / getMilliseconds()) : 0.0f; }

		operator float() const { return m_Time; }
		operator std::string() const { return std::to_string(m_Time); }

	private:
		float m_Time;
	};

	class Timer {
	public:
		explicit Timer(std::string description = "")
			: m_Description(std::move(description)),
			m_Started(false),
			m_Ended(false),
			m_TimeDifference(0.0f)
		{
		}

		~Timer() {
			if (!m_Ended) {
				stop();
			}
		}

		void start() {
			m_Started = true;
			m_Ended = false;
			m_StartTime = std::chrono::high_resolution_clock::now();
		}

		bool isEnded() const { return m_Ended; }

		void stop() {
			if (!m_Started) {
				return;
			}

			auto endTime = std::chrono::high_resolution_clock::now();
			m_TimeDifference = std::chrono::duration<float, std::chrono::milliseconds::period>(endTime - m_StartTime).count();

			m_Ended = true;
			m_Started = false;

			if (!m_Description.empty()) {
				std::cout << m_Description << " took " << m_TimeDifference << " ms" << std::endl;
			}
		}

		float getTimeMilliseconds() const { return m_TimeDifference; }
		float getTimeSeconds() const { return m_TimeDifference / 1000.0f; }

	private:
		bool m_Started;
		bool m_Ended;
		float m_TimeDifference;
		std::string m_Description;
		std::chrono::time_point<std::chrono::high_resolution_clock> m_StartTime;
	};
}