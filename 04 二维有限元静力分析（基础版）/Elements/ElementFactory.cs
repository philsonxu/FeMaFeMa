using System;
using FEM2D.Mesh;

namespace FEM2D.Elements
{
    /// <summary>
    /// 单元工厂：根据 ElementType 创建对应的单元对象。
    /// </summary>
    public static class ElementFactory
    {
        public static FiniteElement Create(ElementType type)
        {
            switch (type)
            {
                case ElementType.CST3: return new CST3Element();
                case ElementType.LT6: return new LT6Element();
                case ElementType.Q4: return new Q4Element();
                default: throw new ArgumentException("Unknown element type: " + type);
            }
        }
    }
}
