#include "Renderer.hpp"

int main()
{
	AR::Renderer renderer{};
	renderer.init(2000, 1400);
	renderer.run();
}