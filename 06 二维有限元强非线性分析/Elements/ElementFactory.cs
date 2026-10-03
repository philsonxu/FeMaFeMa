using System;

namespace NonlinearFEM2D.Elements
{
    public static class ElementFactory
    {
        public static FiniteElement Create(ElementType t)
        {
            switch (t)
            {
                case ElementType.CST3: return new CST3Element();
                case ElementType.LT6: return new LT6Element();
                case ElementType.Q4: return new Q4Element();
                default: return new CST3Element();
            }
        }
    }
}
