//#include "frame_arena.hpp"
//
//namespace AR
//{
//	void* FrameArena::allocate(size_t size, size_t alignment)
//	{
//		size_t current = (size_t)m_DataPtr + m_Offset;
//		size_t aligned = (current + (alignment - 1)) & ~(alignment - 1);
//		size_t nextOffset = (aligned + size) - (size_t)m_DataPtr;
//
//		if (nextOffset > m_ArenaSize)
//		{
//			throw std::bad_alloc();
//		}
//		m_Offset = nextOffset;
//		return (void*)aligned;
//	}
//}