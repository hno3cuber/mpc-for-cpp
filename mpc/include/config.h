#pragma once

constexpr auto N = 20; //预测步长

constexpr auto STANUM = 6; //状态数量
constexpr auto CONNUM = 2; //输入数量

constexpr auto DT = 0.01;
constexpr auto DTT = DT * DT;

constexpr auto COLS = 803;
constexpr auto RAWS = 2;

constexpr auto CONSROWNUM = N * 1; //约束数量