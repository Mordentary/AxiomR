#include "tracy\Tracy.hpp"
// Overload global new
void* operator new(size_t size) {
	void* ptr = std::malloc(size); // Allocate memory.
	if (!ptr) throw std::bad_alloc(); // Handle allocation failure.
	TracyAlloc(ptr, size);           // Notify Tracy.
	return ptr;
}

// Overload global delete
void operator delete(void* ptr) noexcept {
	if (ptr) {
		TracyFree(ptr);             // Notify Tracy.
		std::free(ptr);             // Free memory.
	}
}

// Overload array new
void* operator new[](size_t size) {
	void* ptr = std::malloc(size); // Allocate memory.
	if (!ptr) throw std::bad_alloc(); // Handle allocation failure.
	TracyAlloc(ptr, size);           // Notify Tracy.
	return ptr;
}

// Overload array delete
void operator delete[](void* ptr) noexcept {
	if (ptr) {
		TracyFree(ptr);             // Notify Tracy.
		std::free(ptr);             // Free memory.
	}
}
#include "Renderer.hpp"
int main()
{
	AR::Renderer renderer{};
	renderer.init(1240, 720);
	renderer.run();
}