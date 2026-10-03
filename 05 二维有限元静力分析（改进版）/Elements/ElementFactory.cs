// ElementFactory.cs - 单元工厂：根据 ElementType 枚举创建对应的有限元实例
using System;

namespace FEM2D.Elements
{
    public static class ElementFactory
    {
        public static FiniteElement Create(ElementType type)
        {
            switch (type)
            {
                case ElementType.CST3:
                    return new CST3Element();
                case ElementType.LT6:
                    return new LT6Element();
                case ElementType.Q4:
                    return new Q4Element();
                default:
                    throw new NotSupportedException(
                        string.Format("不支持的单元类型: {0}", type));
            }
        }
    }
}
