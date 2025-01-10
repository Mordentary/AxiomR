#include"thread_pool.hpp"

namespace AR
{
	ThreadPool::ThreadPool(size_t numThreads) : stop(false) {
		for (size_t i = 0; i < numThreads; ++i) {
			workers.emplace_back([this] {
				while (!this->stop.load(std::memory_order_relaxed)) {
					std::function<void()> task;
					{
						std::unique_lock<std::mutex> lock(this->queueMutex);
						this->condition.wait(lock, [this] {
							return this->stop.load(std::memory_order_relaxed) || !this->tasks.empty();
							});
						if (this->stop.load(std::memory_order_relaxed) && this->tasks.empty())
							return;
						task = std::move(this->tasks.front());
						this->tasks.pop();
					}
					task();
				}
				});
		}
	}
	ThreadPool::~ThreadPool() {
		stop.store(true, std::memory_order_relaxed);
		condition.notify_all();
		for (std::thread& worker : workers)
			worker.join();
	}
}