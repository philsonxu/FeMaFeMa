using System;

namespace MultiPhysicsFEM2D.Elements
{
    public static class ElementFactory
    {
        public static FiniteElement Create(ElementType type)
        {
            switch (type)
            {
                case ElementType.CST3: return new CST3Element();
                case ElementType.LT6: return new LT6Element();
                case ElementType.Q4: return new Q4Element();
                default: throw new NotSupportedException("Unsupported element type: " + type);
            }
        }
    }
}
