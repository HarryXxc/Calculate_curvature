# Calculate_curvature
轮廓曲线的曲率计算
代码中有两种方案
1. 如果调用第一种函数CTAR_cal()，那么list是一个二维矩阵，第一维是x，序号的坐标[0,1,2,3,4...n];第二维是曲线的纵坐标点
2. 如果调用第二种函数cal_RTPDA(n, x, y)，n就是你设置的领域大小，不过x和y是横坐标和纵坐标，大小是2*n+1
