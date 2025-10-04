#include <iostream>
#include <vector>
#include "triangles.h"
#include <set>
#include <algorithm>
#pragma once


//компиляция: ./build/triangles < ./tests/test1_input.txt или же ./build/triangles если ввод с клавиатуры

int main() {
    int N;
    std::cout << "Введите число треугольников: ";
    std::cin >> N;

    std::vector<Triangle3D<double>> triangles(N);
    std::cout << "Введите координаты точек треугольников: " << std::endl;
    for (int i = 0; i < N; i++) {
        std::cin >> triangles[i].v0.x >> triangles[i].v0.y >> triangles[i].v0.z;
        std::cin >> triangles[i].v1.x >> triangles[i].v1.y >> triangles[i].v1.z;
        std::cin >> triangles[i].v2.x >> triangles[i].v2.y >> triangles[i].v2.z;
    }

    bool found = false;
    std::set<int> triangles_set;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            auto res = TriangleIntersection(triangles[i], triangles[j]);
            if (res.hit) {
                triangles_set.insert(i);
                triangles_set.insert(j);
                //std::cout << i << " " << j << "\n";
                found = true;
            }
        }
    }

    if (!found) {
        std::cout << "Пересечений нет" << std::endl;
    }
    else {
        std::cout << "Итого, вот номера: " << std::endl;
        for (int tri : triangles_set) {
            std::cout << tri << std::endl;
        }
    }


    return 0;
}
