"""
function:计算离散点曲线的曲率
"""
import math
import numpy as np


def CTAR_cal(list):
    """
    利用CTAR的方式计算，具体见论文：重庆大学《几种轮廓曲率估计角点检测算法研究》
    这个方法在计算三角形的外接圆的曲率，效果还算好
    list是个二维矩阵，表示离散点
    """
    qulv = []
    for i in range(5, len(list[:, 0]) - 5, 5):
        x, y = list[i, 0], list[i, 1]  # 当前点的坐标
        x0, y0 = list[i - 5, 0], list[i - 5, 1]
        x1, y1 = list[i + 5, 0], list[i + 5, 1]

        d1 = math.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
        d2 = math.sqrt((x - x0) ** 2 + (y - y0) ** 2)
        d3 = math.sqrt((x - x1) ** 2 + (y - y1) ** 2)

		# 这边直接求三角形的外接圆半径倒数，计算量大
        cos = (d2 ** 2 + d3 ** 2 - d1 ** 2) / (2 * d2 * d3)
        sin = math.sqrt(1 - cos ** 2)
        k = 2 * sin / d1
        # k = d1 / (d2 + d3)   # 论文中用的是这个方法，
        qulv.append(k)

    return qulv


def cal_RTPDA(n, x, y):
    """
    :param n:邻域大小，总的数组长度应该是2*n+1
    :param x: 长度为2*n+1的x轴
    :param y: 长度为2*n+1的y轴
    :return: 曲率大小，和你选的邻域大小有关，具体参考那篇博士论文
    """
    s = [0 for s in range(2 * n + 1)]  # 弧长
    s[n] = 0
    for i in range(n - 1, -1, -1):
        s[i] = s[i + 1] - np.sqrt((x[i] - x[i + 1]) ** 2 + (y[i] - y[i + 1]) ** 2)

    for i in range(n + 1, 2 * n + 1):
        s[i] = s[i - 1] + np.sqrt((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2)

    a11, a12, a22, ex, fx = 0, 0, 0, 0, 0
    b11, b12, b22, ey, fy = 0, 0, 0, 0, 0
    for i in range(0, 2 * n + 1):
        a11 += s[i] ** 2
        a12 += s[i] ** 3
        a22 += s[i] ** 4
        ex += s[i] * (x[i] - x[n])
        fx += (s[i] ** 2) * (x[i] - x[n])

        b11 += s[i] ** 2
        b12 += s[i] ** 3
        b22 += s[i] ** 4
        ey += s[i] * (y[i] - y[n])
        fy += (s[i] ** 2) * (y[i] - y[n])

    dx = a11 * a22 - a12 ** 2
    af0 = x[n]
    af1 = (a22 * ex - a12 * fx) / dx
    af2 = (a11 * fx - a12 * ex) / dx

    dy = b11 * b22 - b12 ** 2
    bf0 = y[n]
    bf1 = (b22 * ey - b12 * fy) / dy
    bf2 = (b11 * fy - b12 * ey) / dy

    RTPDA = 0
    for i in range(0, 2 * n + 1):
        if i == n:  # 这里很重要，一定要排除这个点
            continue
        temp1 = bf1 * (x[i] - x[n]) - af1 * (y[i] - y[n])
        temp2 = np.sqrt((x[i] - x[n]) ** 2 + (y[i] - y[n]) ** 2)
        RTPDA += abs(temp1 / temp2)

    RTPDA = RTPDA * (1 / np.sqrt(af1 ** 2 + bf1 ** 2))

    return RTPDA
