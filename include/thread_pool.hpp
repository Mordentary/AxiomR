#pragma once
#include <queue>
#include <memory>
#include <mutex>
#include <future>
#include <type_traits>
namespace AR
{
	class ThreadPool {
	public:
		ThreadPool(size_t numThreads);
		~ThreadPool();
		template <class F, class... Args>
		auto enqueue(F&& f, Args&&... args) -> std::future<typename std::invoke_result<F, Args...>::type>
		{
			using return_type = typename std::invoke_result<F, Args...>::type;

			auto task = std::make_shared< std::packaged_task<return_type()> >(std::bind(std::forward<F>(f), std::forward<Args>(args)...));

			std::future<return_type> res = task->get_future();
			{
				std::unique_lock<std::mutex> lock(queueMutex);
				if (stop.load(std::memory_order_relaxed))
					throw std::runtime_error("enqueue on stopped ThreadPool");
				tasks.emplace([task]() { (*task)(); });
			}
			condition.notify_one();
			return res;
		}
	private:
		std::vector<std::thread> workers;
		std::queue<std::function<void()>> tasks;
		std::mutex queueMutex;
		std::condition_variable condition;
		std::atomic<bool> stop;
	};
}