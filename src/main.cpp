#include "Renderer.hpp"

int main()
{
	AR::Renderer renderer{};
	renderer.init(1000, 1000);
	renderer.run();
}