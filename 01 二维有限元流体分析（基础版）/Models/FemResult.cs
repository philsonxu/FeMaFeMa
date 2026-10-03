using System;

namespace Fem2DFluid.Models
{
    /// <summary>
    /// 求解结果摘要
    /// </summary>
    public class FemResult
    {
        public bool Success;
        public string Message;
        public int NodeCount;
        public int ElementCount;
        public double PhiMin;
        public double PhiMax;
        public double VmagMin;
        public double VmagMax;
        public TimeSpan SolveTime;
    }
}
