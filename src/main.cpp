#include "tracy\Tracy.hpp"
#include <cstdlib>
#include <new>

// Overload global new
void* operator new(size_t size) {
	void* ptr = std::malloc(size);     
	if (!ptr) throw std::bad_alloc(); 
	TracyAlloc(ptr, size);           
	return ptr;
}

// Overload global delete
void operator delete(void* ptr) noexcept {
	if (ptr) {
		TracyFree(ptr);             
		std::free(ptr);            
	}
}

// Overload array new
void* operator new[](size_t size) {
	void* ptr = std::malloc(size);    
	if (!ptr) throw std::bad_alloc(); 
	TracyAlloc(ptr, size);           
	return ptr;
}

// Overload array delete
void operator delete[](void* ptr) noexcept {
	if (ptr) {
		TracyFree(ptr);             
		std::free(ptr);            
	}
}
#include "Renderer.hpp"
int main()
{
	AR::Renderer renderer{};
	renderer.init(1240, 720);
	renderer.run();
}