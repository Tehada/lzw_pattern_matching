# lzw pattern matching

## Установка зависимостей

Данный проект использует пакетный менеджер [vcpkg](https://vcpkg.io/en/). Нужно установить следующие зависимости, запустив у себя в терминале команды:
```
vcpkg install spdlog
vcpkg install gtest
```

## Запуск

Для сборки проекта используется cmake. Следующие команды скомпилируют код в исполняемый файл:
```
cmake -B build -S . "-DCMAKE_TOOLCHAIN_FILE=/opt/vcpkg/scripts/buildsystems/vcpkg.cmake" -DCMAKE_CXX_COMPILER=clang++
cmake --build build
```

Опцию `DCMAKE_TOOLCHAIN_FILE`, возможно, придется выставить под своё окружение. Опцию `DCMAKE_CXX_COMPILER` можно не выставлять, если в окружении нет компилятора `clang`.

После запуска этих команд нужно подготовить файл со сжатым текстом:
```
clang++ -std=c++20 lzw.cpp && ./a.out -c input2 output2
```

Команда для запуска поиска по сжатому тексту:
```
./build/main output2
```
