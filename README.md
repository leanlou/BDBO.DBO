# luck_dog
# xdm,我是中国人，我就不写英文了哈

# 蜣螂优化（Dung Beetle Optimization，DBO）算法是一种适用于求解高维非凸优化问题的启发式算法，但其随机搜索能力存在缺陷并且而容易陷入局部最优。为了增强其搜索能力，同时提升算法的优化精度，研究选择拓宽种群探索范围，引入位置更新辅助策略，同时增强其局部开发能力。

# 滚球蜣螂种群搜索范围依据控制因子R呈线性缩减，研究替换其更新公式为抛物线形变化来扩展算法搜索范围；为帮助算法逃离局部最优，研究引入Levy飞行机制辅助滚球蜣螂进行远距离位置更新，同时增强其全局范围探索能力；育雏蜣螂与觅食蜣螂的权重因子β随机性过高不利于局部寻优，研究控制随机权重β呈高斯分布来提升算法的优化精度。

# PS：高斯分布权重因子（GDBO）造成的随机扰动对DBO算法性能提升具有显著作用，Levy飞行与控制因子R替换策略并无显著效果。
# GDBO的性能其实优于DBO与BDBO
# g设定均值为0，方差为1，效果不错，继续调参需要做田口实验（懒得做）

# function （Multi_Convergence_curve, Multi_time, bestfit） = bdbo(Z,Popsize, dim, ub, lb, fobj, Max_iteration, Multi_Num,func_num)
# 其中， Multi_Convergence_curve为收敛曲线，Multi_time为运行时间，bestfit为算法找到的最优值
# Z为随机种子数，Popsize为种群数, dim为维度, ub上界upper, lb下界lower, fobj函数定义，main主函数运行中前面要加@,
#  Max_iteration最大迭代次数, Multi_Num最大运行次数，次数设置为1,func_num为CEC函数序号。


# 注意代码中局部最优解的编码方式，会对算法性能产生较大影响（亲测）
