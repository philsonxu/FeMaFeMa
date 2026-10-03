using System;

namespace Fem2DFluid.Models
{
    /// <summary>
    /// 材料参数：二维势流/导热问题中，相当于渗透率k或导热系数λ
    /// 控制方程：∇·(k ∇φ) = f
    /// 势流中：φ为速度势，k=1，速度 u=∂φ/∂x, v=∂φ/∂y
    /// 渗流中：φ为水头，k为渗透率
    /// 导热中：φ为温度，k为导热系数
    /// </summary>
    public class Material
    {
        public int Id;
        public string Name;
        public double K; // 渗透率/导热系数（常数各向同性）
        public double Source; // 源项f（单位体积产热/单位体积流出/源汇）

        public Material()
        {
        }

        public Material(int id, string name, double k)
        {
            Id = id;
            Name = name;
            K = k;
            Source = 0.0;
        }
    }
}
