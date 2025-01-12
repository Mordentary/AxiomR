//#pragma once
//#include<memory>
//namespace AR
//{
//	class FrameArena
//	{
//	public:
//		FrameArena(size_t sizeInBytes) : m_ArenaSize(sizeInBytes), m_Offset(0), m_DataPtr(::operator new(sizeInBytes););
//		~FrameArena() { ::operator delete(m_DataPtr); };
//		void* allocate(size_t size, size_t alignment = alignof(std::max_align_t));
//		inline void reset() { m_Offset = 0; }
//	private:
//		size_t m_Offset = 0;
//		size_t m_ArenaSize = 0;
//		void* m_DataPtr = nullptr;
//	};
//}