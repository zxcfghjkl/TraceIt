// Для загрузки модели, избавиться.
#include <strstream>
// Для дебаггинга
#include <iostream>
// Работа с файлами
#include <fstream>
// Триганометрия, степени и корни
#include <cmath>
// Вектор
#include <vector>
// Строка (тип данных)
#include <string>

#ifdef _WIN32 // (Windows)
// Аналог usleep для Windows
	#include <windows.h>
	#define usleep(ms) Sleep((ms) / 1000) // В Windows миллисекунды
#else // (MacOS / Linux / BSD / Android, и т. п.)
// Для usleep()
	#include <unistd.h>
#endif
#ifdef _WIN32
// Вывод графики
#else
    #include <ncurses.h> // Для NCurses
#endif

using namespace std;

int height, width;

// Вектор с 3 элементами
struct vec3d
{
    float x, y, z;
};

// 4x4 матрица
struct matrx4d
{
    float m[4][4]; // Двумерный массив для удобства работы
};

//enum vertDir = {clockwise, antiClockwise};
// Треугольник
struct triangle
{
    vec3d p[3];
    // vertDir dir = clockwise;
};

// Сетка
struct mesh {
    vector<triangle> triangles; // Массив треугольников, составляющих сетку
    vec3d boxSize; // Размер модели (и хитбокса)
};

// Коллайдер
struct collider {
    vec3d center;
    vec3d size; // Для AABB (коробки)
    float radius; // Для сферы
    bool isSphere;
    bool isStatic;
    float rest;  // Упругость
};

// Структура, описавающая физическое тело
struct object
{
    vec3d pos;  // Позиция обьекта
    vec3d vel; // Скорость
    vec3d force; // Сила
    float mass; // Масса
    mesh model; // Сетка
    matrx4d transform; // Матрица трансформации обьекта, добавить rotMatrix
    float Cd; // Сопротивление формы
    float area;  // Площадь поперечного сечения
    bool stationary; // Не подвергается физике
    vector<vec3d> trajectory; // История позиций
    int maxTrajectoryPoints = 300; // Макс. точек траектории
    bool trajWritten = false;
    collider col;
};


// , Камера
struct camera
{
    object physics;  // Физика камеры
    vec3d position;  // Позиция камеры
    vec3d direction; // Направление взгляда
    vec3d up;        // Вектор "вверх" камеры
};

// Совокупность всех физических тел в мире и параметров мира
struct world
{
    vector<object> objects; // Массив тел
    float globalScale;      // Коэффициент размера
    float atmDensity;        // Плотность атмосферы
    matrx4d projMatrix;     // Матрица проекции
    matrx4d viewMatrix;     // Матрица вида
    float g;                // УСП
};

// Кнопка в меню
enum UIElementType { BUTTON, CHECKBOX };

struct button {
    UIElementType type; // Добавляем тип элемента
    int x, y;
    union {
        float* targetVar;   // Для обычных кнопок
        bool* targetBool;  // Для чекбоксов / флажков
    };
    std::string label;
    int minVal;
    int maxVal;
    float step;
};

// Начальная настройка камеры
camera cam = { {},
    {0.0f, 0.0f, -20.0f},  // Позиция камеры
    {0.0f, 0.0f, 1.0f},   // Направление взгляда (смотрит по оси Z)
    {0.0f, -1.0f, 0.0f}    // Вектор "вверх" (по оси Y)
};

// Функция, использующая алгоритм Брезенхема для рисования линий
void drawLine(int x1, int y1, int x2, int y2) {
    // Проверка выхода за пределы экрана
    //if ((x1 < 0 && x2 < 0) || (x1 >= width && x2 >= width) ||
        //(y1 < 0 && y2 < 0) || (y1 >= height && y2 >= height)) {
        //return;
    //}

    // Сам Брезенхем
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int sx = (x1 < x2) ? 1 : -1;
    int sy = (y1 < y2) ? 1 : -1;
    int err = dx - dy;

    while (true) {
        // Проверка границ перед рисованием
        if (x1 >= 0 && x1 < width && y1 >= 0 && y1 < height) {
            mvaddch(y1, x1, '*');
        }

        if (x1 == x2 && y1 == y2) break;

        int e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x1 += sx;
        }
        if (e2 < dx) {
            err += dx;
            y1 += sy;
        }
    }
}

//  OneLoneCoder'a
mesh loadFromObjectFile(string sFilename)
	{
	    mesh outputMesh;

        // Открываем файл
		ifstream f(sFilename);

		// Локальный кеш вершин
		vector<vec3d> verts;

        // Пока не закончился
		while (!f.eof())
		{
			char line[128];
			f.getline(line, 128);

			strstream s;
			s << line;

			char junk;

            // Если 1 символ "v" - вершина
			if (line[0] == 'v')
			{
				vec3d v;
				s >> junk >> v.x >> v.y >> v.z;
				verts.push_back(v);
			}

            // Это -
			if (line[0] == 'f')
			{
				int f[3];
				s >> junk >> f[0] >> f[1] >> f[2];
				outputMesh.triangles.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
			}
		}

		return outputMesh;
	}

    // Произведения векторов
// Скалярное
float dot(const vec3d &a, const vec3d &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Векторное
vec3d cross(const vec3d &a, const vec3d &b)
{
    vec3d result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

// И нормализация
vec3d normalize(const vec3d &v)
{
    float length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (length != 0.0f)
    {
        vec3d result;
        result.x = v.x / length;
        result.y = v.y / length;
        result.z = v.z / length;
        return result;
    }
    return v;  // Возвращаем исходный вектор, если его длина == 0
}

// Вычитание
vec3d vecSubtract(const vec3d &a, const vec3d &b)
{
    vec3d result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return result;
}

// Сложение
vec3d vecAdd(const vec3d &a, const vec3d &b)
{
    vec3d result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}

// Умножение и деление вектора на скаляр
vec3d multVecScal(const vec3d &vec, const float flt)
{
    vec3d result;
    result.x = vec.x * flt;
    result.y = vec.y * flt;
    result.z = vec.z * flt;
    return result;
}

vec3d divVecScal(const vec3d &vec, const float flt)
{
    vec3d result;
    result.x = vec.x / flt;
    result.y = vec.y / flt;
    result.z = vec.z / flt;
    return result;
}

// Умножение двух матриц 4x4
matrx4d matMult(const matrx4d& a, const matrx4d& b) {
    matrx4d result = {0};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                result.m[i][j] += a.m[i][k] * b.m[k][j];
            }
        }
    }
    return result;
}

// Функция для вычисления нормали треугольника
vec3d calculateNormal(const triangle &tri)
{
    vec3d edge1 = vecSubtract(tri.p[2], tri.p[1]);
    vec3d edge2 = vecSubtract(tri.p[1], tri.p[0]);
    return normalize(cross(edge1, edge2));
}

// Бек-фейс куллинг, нужно clockwise, !clockwise для треугольника
bool isTriangleVisible(const triangle &tri, const camera& c)
{
    vec3d normal = calculateNormal(tri);
    vec3d viewDirection = vecSubtract(c.position, tri.p[0]);

    // Нормализуем векторы
    viewDirection = normalize(viewDirection);

    // Треугольник видим, если нормаль направлена ПРОТИВ вектора взгляда
    return dot(normal, viewDirection) < 0.0f;
}

// Функция умножения матрицы на вектор
vec3d matVecMult(const matrx4d &mat, const vec3d &vec)
{
    float w = mat.m[3][0] * vec.x + mat.m[3][1] * vec.y + mat.m[3][2] * vec.z + mat.m[3][3];
    if (w == 0.0f) w = 1.0f; // Предотвращаем деление на ноль
    vec3d result;
    result.x = (mat.m[0][0] * vec.x + mat.m[0][1] * vec.y + mat.m[0][2] * vec.z + mat.m[0][3]) / w;
    result.y = (mat.m[1][0] * vec.x + mat.m[1][1] * vec.y + mat.m[1][2] * vec.z + mat.m[1][3]) / w;
    result.z = (mat.m[2][0] * vec.x + mat.m[2][1] * vec.y + mat.m[2][2] * vec.z + mat.m[2][3]) / w;
    return result;
}

// -1 - 1 >> 0 - конец экрана
vec3d scaleToScreen(const vec3d &point)
{
    vec3d scaled;
    scaled.x = (point.x + 1.0f) * 0.5f * (width - 1);  // Преобразование [-1, 1] в [0, width]
    scaled.y = (1.0f - point.y) * 0.5f * (height - 1); // Инверсия Y и преобразование [-1, 1] в [0, height]
    scaled.z = point.z; // Z оставляем как есть
    return scaled;
}

// Функция для рисования треугольника с учетом всех преобразований
void drawTriangle(triangle tri, const matrx4d& projMatrix, const matrx4d& viewMatrix,
                 const matrx4d& scaleMatrix, const matrx4d& rotMatrix, const matrx4d& transMatrix) {
    bool allOutside = true;
    for (int i = 0; i < 3; i++) {
        // Масштаб > Поворот > Перемещение > Вид > Проекция
        vec3d point = tri.p[i];
        point = matVecMult(scaleMatrix, point);
        point = matVecMult(rotMatrix, point);
    	point = matVecMult(transMatrix, point);
        point = matVecMult(viewMatrix, point);
        point = matVecMult(projMatrix, point);
        tri.p[i] = scaleToScreen(point);
        // Проверка попадания в экран
        if (tri.p[i].x >= 0 && tri.p[i].x < width &&
            tri.p[i].y >= 0 && tri.p[i].y < height) {
            allOutside = false;
        }
    }

    if(allOutside) return;

    drawLine(tri.p[0].x, tri.p[0].y, tri.p[1].x, tri.p[1].y);
    drawLine(tri.p[1].x, tri.p[1].y, tri.p[2].x, tri.p[2].y);
    drawLine(tri.p[2].x, tri.p[2].y, tri.p[0].x, tri.p[0].y);

    //move(10, 0);

    //printw("Point 1: (x%.2f, y%.2f) \n", tri.p[0].x, tri.p[0].y);
    //printw("Point 2: (x%.2f, y%.2f) \n", tri.p[1].x, tri.p[1].y);
    //printw("Point 3: (x%.2f, y%.2f) \n", tri.p[2].x, tri.p[2].y);
}

// Функция для масштабирования
triangle scaleTri(const triangle &tri, float scale)
{
    triangle result;
    for (int i = 0; i < 3; i++)
    {
        // Масштабируем все координаты
        result.p[i].x = tri.p[i].x * scale;
        result.p[i].y = tri.p[i].y * scale;
        result.p[i].z = tri.p[i].z * scale; // Масштабируем Z тоже
    }
    return result;
}

// Функция для создания матрицы масштаба
matrx4d createScaleMatrix(float scaleX, float scaleY, float scaleZ)
{
    matrx4d scaleMat = {{{scaleX, 0.0f, 0.0f, 0.0f},
                        {0.0f, scaleY, 0.0f, 0.0f},
                        {0.0f, 0.0f, scaleZ, 0.0f},
                        {0.0f, 0.0f, 0.0f, 1.0f}}};
    return scaleMat;
}

// Функция для создания матрицы вида
matrx4d createViewMatrix(const camera& cam)
{
    // Векторы для камеры
    vec3d zAxis = normalize(cam.direction);
    vec3d xAxis = normalize(cross(cam.up, zAxis));
    vec3d yAxis = cross(zAxis, xAxis);

    // Матрица вида
    matrx4d viewMatrix = {{
    {xAxis.x, xAxis.y, xAxis.z, -dot(xAxis, cam.position)},
    {yAxis.x, yAxis.y, yAxis.z, -dot(yAxis, cam.position)},
    {zAxis.x, zAxis.y, zAxis.z, -dot(zAxis, cam.position)},
    {0.0f, 0.0f, 0.0f, 1.0f}
}};

    return viewMatrix;
}

// Функции для создания матрицы поворота

matrx4d createRotationX(float angle) {
    float c = cos(angle), s = sin(angle);
    return {{
        {1, 0,  0, 0},
        {0, c, -s, 0},
        {0, s,  c, 0},
        {0, 0,  0, 1}
    }};
}

matrx4d createRotationY(float angle) {
    float c = cos(angle), s = sin(angle);
    return {{
        {c,  0, s, 0},
        {0,  1, 0, 0},
        {-s, 0, c, 0},
        {0,  0, 0, 1}
    }};
}

matrx4d createRotationZ(float angle) {
    float c = cos(angle), s = sin(angle);
    return {{
        {c, -s, 0, 0},
        {s,  c, 0, 0},
        {0,  0, 1, 0},
        {0,  0, 0, 1}
    }};
}

// Функция для создания матрицы трансляции (перемещения)
matrx4d createTranslationMatrix(const float& transX, const float& transY, const float& transZ)
{
        matrx4d transMatrix = {{
                {1.0f, 0.0f, 0.0f, transX},
                {0.0f, 1.0f, 0.0f, transY},
                {0.0f, 0.0f, 1.0f, transZ},
                {0.0f, 0.0f, 0.0f, 1.0f}
        }};

        return transMatrix;
}

// Функция для создания матрицы трансляции проекции
matrx4d createProjMatrix(const float& aspect, const float& f)
{
        const float near = 0.1f;
        const float far = 100.0f;
        const float fov = f * M_PI / 180.0f;
        const float fovRad = 1.0f / tan(fov * 0.5f);
        matrx4d projMatrix = {{ // 2.519 - соотношение сторон символа
          {fovRad * aspect * 2.519f, 0.0f, 0.0f, 0.0f},
          {0.0f,            fovRad,    0.0f,                 0.0f},
          {0.0f,            0.0f,      (far + near)/(far - near), 1.0f},
          {0.0f,            0.0f,      (-2.0f * far * near)/(far - near), 0.0f}
          }};

        return projMatrix;
}

// Клиппинг
mesh clipTri(const triangle &tri)
{
    mesh clipped;
    if(abs(tri.p[0].y) - height > 0 && abs(tri.p[0].x) - width > 0 && abs(tri.p[1].y) - height > 0 && abs(tri.p[1].x) - width > 0 && abs(tri.p[2].y) - height > 0 && abs(tri.p[2].x) - width > 0) return clipped = {{tri}};
    //else if()
}

// Рисуем сетку
void drawMesh(const mesh& inputMesh, const matrx4d& projMatrix, const matrx4d& viewMatrix, float globalScale,
 		const matrx4d& rotMatrix, const matrx4d& transMatrix, const camera& c) {
    // Создаём матрицу масштаба
    matrx4d scaleMatrix = createScaleMatrix(globalScale, globalScale, globalScale);

    // Проходимся по всем треугольникам в сетке
    for (const triangle& trian : inputMesh.triangles) {
        //float y = trian.p[2].y - trian.p[1].y;
        //float step = abs(y - trian.p[0].y) / 10;
        //for(int i = 0; i < 10; i++)
        //{
           // while(abs(y) < trian.p[0].y)
            //{
            //drawLine(int x1, int y1, int x2, int y2);
            //y+= step;
           // }
       //}

        // Проверяем видимость треугольника
        //if (isTriangleVisible(trian, c)) {
            drawTriangle(trian, projMatrix, viewMatrix, scaleMatrix, rotMatrix, transMatrix); // Отрисовка видимого треугольника
        //}
    }
}

// Аналог cin >> из стандартной библиотеки
string consoleIn() {
    string input;
    // Проверка if(), else() return "Error"
    keypad(stdscr, TRUE);

    int ch;
    while ((ch = getch()) != '\n' && ch != KEY_ENTER) {
        if (ch == KEY_BACKSPACE || ch == 127 || ch == 8) { // Обработка Backspace
            if (!input.empty()) {
                input.pop_back();
                // Удалить символ с экрана
                addstr("\b \b");
                refresh();
            }
        } else if (isprint(ch)) { // Проверка на печатный символ
            input.push_back(ch);
            addch(ch); // Отобразить символ на экране
            refresh();
        }
    }
    return input;
}

// Создание паралелепипеда
mesh createRectMesh(const float &sX, const float &sY, const float &sZ) {
    mesh rectMesh;

    // Вершины
    vec3d vertices[8] = {
    {-0.5f * sX, -0.5f * sY, -0.5f * sZ}, // 0
    {0.5f * sX, -0.5f * sY, -0.5f * sZ}, // 1
    {0.5f * sX, 0.5f * sY, -0.5f * sZ}, // 2
    {-0.5f * sX, 0.5f * sY, -0.5f * sZ}, // 3
    {-0.5f * sX, -0.5f * sY, 0.5f * sZ}, // 4
    {0.5f * sX, -0.5f * sY, 0.5f * sZ}, // 5
    {0.5f * sX, 0.5f * sY, 0.5f * sZ}, // 6
    {-0.5f * sX, 0.5f * sY, 0.5f * sZ}  // 7
    };

    // Определяем треугольники для каждой грани
    triangle cubeTriangles[12] = {
        // Задняя грань
        {vertices[0], vertices[1], vertices[2]},
        {vertices[2], vertices[3], vertices[0]},

        // Передняя грань
        {vertices[4], vertices[5], vertices[6]},
        {vertices[6], vertices[7], vertices[4]},

        // Левая грань
        {vertices[0], vertices[3], vertices[7]},
        {vertices[7], vertices[4], vertices[0]},

        // Правая грань
        {vertices[1], vertices[5], vertices[6]},
        {vertices[6], vertices[2], vertices[1]},

        // Нижняя грань
        {vertices[0], vertices[1], vertices[5]},
        {vertices[5], vertices[4], vertices[0]},

        // Верхняя грань
        {vertices[2], vertices[6], vertices[7]},
        {vertices[7], vertices[3], vertices[2]}
    };

    // Добавляем треугольники в mesh
    for (const auto& tri : cubeTriangles) { // Для каждого треугольника в эррее "cubeTriangles"
        rectMesh.triangles.push_back(tri);
    }

    return rectMesh;
}

// Нормаль коробки
vec3d calculateAABBNormal(const collider& a, const collider& b) {
    vec3d delta = vecSubtract(b.center, a.center);
    vec3d overlap = {
        (a.size.x + b.size.x) - abs(delta.x),
        (a.size.y + b.size.y) - abs(delta.y),
        (a.size.z + b.size.z) - abs(delta.z)
    };

    // Находим ось с минимальным перекрытием
    if(overlap.x < overlap.y && overlap.x < overlap.z) {
        return {delta.x > 0 ? 1.0f : -1.0f, 0, 0};
    } else if(overlap.y < overlap.z) {
        return {0, delta.y > 0 ? 1.0f : -1.0f, 0};
    } else {
        return {0, 0, delta.z > 0 ? 1.0f : -1.0f};
    }
}

// Обновление матрицы трансляции в обьекте
void updateObjectTransform(object& obj) {
    obj.transform = createTranslationMatrix(obj.pos.x, obj.pos.y, obj.pos.z);
}


// Функция проверки коллизий Коробка-Коробка
bool checkAABBCollision(const collider& a, const collider& b) {
    return (abs(a.center.x - b.center.x) <= (a.size.x + b.size.x) && // Расстояние между центрами a и b меньше/равно сумме их размеров
            abs(a.center.y - b.center.y) <= (a.size.y + b.size.y) &&
            abs(a.center.z - b.center.z) <= (a.size.z + b.size.z));
}

// Функция проверки коллизий Сфера-Сфера
bool checkSphereCollision(const collider& a, const collider& b) {
    vec3d delta = vecSubtract(a.center, b.center); // ~ Тоже самое
    float distanceSq = dot(delta, delta);
    float radiusSum = a.radius + b.radius;
    return distanceSq <= (radiusSum * radiusSum);
}

// Обновляем БАК
void updateCollider(object& obj) {
    obj.col.center = obj.pos;

    // Если не сфера
    //if(!obj.col.isSphere) {
        //obj.col.size = multVecScal(obj.model.boxSize, obj.transform.m[0][0]);
    //}
}

// Коллизии, коллизии, коллизии...
void resolveCollisions(world& w) {
    for(size_t i = 0; i < w.objects.size(); ++i) {
        for(size_t j = i+1; j < w.objects.size(); ++j) {
            object& a = w.objects[i];
            object& b = w.objects[j];

            if(a.col.isStatic && b.col.isStatic) continue;
            if(a.mass <= 0 || b.mass <= 0) continue;

            bool collision = false;
            vec3d normal;

            if(a.col.isSphere && b.col.isSphere) {
                collision = checkSphereCollision(a.col, b.col);
                if(collision) normal = normalize(vecSubtract(b.col.center, a.col.center));
            } else {
                collision = checkAABBCollision(a.col, b.col);
                if(collision) normal = calculateAABBNormal(a.col, b.col);
            }

            if(collision) {
                // Относительная скорость
                vec3d relVelocity = vecSubtract(b.vel, a.vel);
                float velAlongNormal = dot(relVelocity, normal);

                if(velAlongNormal > 0) continue;

                // Импульс
                float e = std::min(a.col.rest, b.col.rest);
                float j = -(1 + e) * velAlongNormal;
                j /= (1/a.mass + 1/b.mass);

                // Применение импульсов
                vec3d impulse = multVecScal(normal, j);

                if(!a.col.isStatic) {
                    a.vel = vecSubtract(a.vel, multVecScal(impulse, 1/a.mass));
                }
                if(!b.col.isStatic) {
                    b.vel = vecAdd(b.vel, multVecScal(impulse, 1/b.mass));
                }

                // Коррекция позиции (менее 1% от перекрытия)
                const float penetrationSlop = 0.01f;
                const float percent = 0.2f;
                vec3d correction = multVecScal(normal, percent * penetrationSlop);

                if(!a.col.isStatic) {
                    a.pos = vecSubtract(a.pos, multVecScal(correction, 1/a.mass));
                }
                if(!b.col.isStatic) {
                    b.pos = vecAdd(b.pos, multVecScal(correction, 1/b.mass));
                }
            }
        }
    }
}

// Обновление физики в мире
void updatePhysics(world& w, ofstream& coords) {
    const vec3d gravity = {0, -9.8, 0}; // Гравитация должна быть отрицательной по Y
    float dt = 0.01666f;

    // Открываем файл один раз для всех объектов
    static bool headerWritten = false;
    coords.open("Координаты.txt", ios::app); // Режим добавления

    // Записываем заголовок только один раз
    if (!headerWritten) {
        coords << "Координаты в формате (x, y, z, время):\n";
        headerWritten = true;
    }

    for(auto& obj : w.objects) {
        if(!obj.stationary)
        {
            // Рассчитываем скорость^2
            float speedSq = dot(obj.vel, obj.vel);

            // Нормализуем скорость для направления
            vec3d velDir = {0, 0, 0};
            if(speedSq > 0.0001f) { // Избегаем деления на ноль
                float invSpeed = 1.0f / sqrt(speedSq); // Скорость, противоположная по длине
                velDir = {obj.vel.x * invSpeed, obj.vel.y * invSpeed, obj.vel.z * invSpeed};
            }

            // Формула силы сопротивления: Fd = -0.5 * p * v^2 * Cd * A * направление_скорости
            vec3d Fd = multVecScal(velDir, -0.5f * w.atmDensity * obj.Cd * obj.area * speedSq);

            // Суммируем все силы (Fd + Fgravity + Fother)
            vec3d totalForce = vecAdd(Fd, multVecScal(gravity, obj.mass));
            totalForce = vecAdd(totalForce, obj.force); // Добавляем другие силы

            // Интегрируем
            vec3d acceleration = divVecScal(totalForce, obj.mass);
            obj.vel = vecAdd(obj.vel, multVecScal(acceleration, dt));
            obj.pos = vecAdd(obj.pos, multVecScal(obj.vel, dt));

            // Обновляем коллайдеры
            updateCollider(obj);

            // Добавляем позицию в историю
            obj.trajectory.push_back(obj.pos);

            // Записываем траекторию, если объект ещё не записан
            if(obj.trajectory.size() == obj.maxTrajectoryPoints && !obj.trajWritten) {
                    coords << "\nТраектория объекта [" << "obj.id" << "]:\n"; // Добавляем ID объекта
                for(const auto& point : obj.trajectory) {
                    coords << point.x << "(м)\t" << point.y << "(м)\t" << point.z << "(м)\n";
                }
                obj.trajWritten = true;
            }

            coords.close(); // Закрываем файл после обработки всех объектов

            // Ограничиваем длину траектории
            if(obj.trajectory.size() > obj.maxTrajectoryPoints) {
                obj.trajectory.erase(obj.trajectory.begin());
            }
            // Обновляем матрицу
            updateObjectTransform(obj);
            obj.force = {0, 0, 0}; // Сброс внешних сил
        }
    }
    // Проверяем и разрешаем коллизии
    resolveCollisions(w);
}

// Рисование траектории тела
void drawTrajectory(object& obj,
                   const matrx4d& projMatrix,
                   const matrx4d& viewMatrix,
                   float globalScale) { //const int& frames
    if(obj.trajectory.size() < 2) return;

    matrx4d scaleMatrix = createScaleMatrix(globalScale, globalScale, globalScale);

    // Преобразуем все точки в экранные координаты
    vector<vec3d> screenPoints;
    for(const auto& point : obj.trajectory) {
        vec3d transformed = matVecMult(scaleMatrix, point);
        transformed = matVecMult(viewMatrix, transformed);
        transformed = matVecMult(projMatrix, transformed);
        screenPoints.push_back(scaleToScreen(transformed));
        }

    // Рисуем линии между точками
    for(size_t i = 1; i < screenPoints.size(); i++) {
        vec3d prev = screenPoints[i-1];
        vec3d curr = screenPoints[i];
        drawLine(prev.x, prev.y, curr.x, curr.y);
    }
}

// Рисование всех обьектов и их траектории в мире
void drawWorld(world &w, matrx4d &rotMatrix, const camera& cam)
{
    for (auto &b : w.objects) { // Для каждого тела в массиве "theWorld.objects"
        // Рисование обьекта
        // matrx4d tfMatrix = matMult(scaleMatrix, (matMult(rotMatrix, (matMult(obj.transform, (matMult(w.projMatrix, w.viewMatrix))))))
        drawMesh(b.model, w.projMatrix, w.viewMatrix,
                 w.globalScale, rotMatrix, b.transform, cam);

        // Отрисовка траектории
        drawTrajectory(b, w.projMatrix,
                       w.viewMatrix, w.globalScale);
    }
}

// Создание окна
WINDOW* createWindow(vector<button>& buttons, int winH, int winW,
                      int startY = -1, int startX = -1,
                      int textClr = COLOR_WHITE, int bgClr = COLOR_BLACK,
                      const char* title = nullptr, bool border = true) {
    // Проверка цветовой поддержки
    if (!has_colors()) {
        throw runtime_error("Terminal does not support colors");
    }

    // Создаем цветовые пары
    init_pair(1, textClr, bgClr);
    init_pair(2, COLOR_BLACK, COLOR_WHITE);

    // Позиционирование окна
    if (startY < 0) startY = (LINES - winH) / 2;
    if (startX < 0) startX = (COLS - winW) / 2;

    WINDOW* win = newwin(winH, winW, startY, startX);
    wbkgd(win, COLOR_PAIR(1));
    wattron(win, COLOR_PAIR(1));

    if (border) {
        box(win, ACS_VLINE, ACS_HLINE);
        if (title) {
            mvwprintw(win, 0, 2, "[ %s ]", title);
        }
    }

    for (auto& btn : buttons) {
    if(btn.type == CHECKBOX) {
        mvwprintw(win, btn.y, btn.x, "[%c] %s",
                *btn.targetBool ? 'X' : ' ',
                btn.label.c_str());
    }
    else {
        mvwprintw(win, btn.y, btn.x, "< %s: %.2f >",
                btn.label.c_str(), *btn.targetVar);
    }
}

    keypad(win, TRUE);
    wrefresh(win);
    return win;
}

// Обновление кнопки
void updateButton(WINDOW* win, button& btn, bool selected) {
    if (selected) wattron(win, COLOR_PAIR(2));

    if(btn.type == CHECKBOX) {
        // Отрисовка чекбокса: [X] или [ ]
        mvwprintw(win, btn.y, btn.x, "[%c] %s",
                *btn.targetBool ? 'X' : ' ',
                btn.label.c_str());
    }
    else {
        // Отрисовка обычной кнопки
        mvwprintw(win, btn.y, btn.x, "< %s: %.2f >",
                btn.label.c_str(), *btn.targetVar);
    }

    if (selected) wattroff(win, COLOR_PAIR(2));
    wrefresh(win);
}

// Управление вводом, switch() - не лучшее решение
void handleInput(int &ch, camera &cam, world &w,
float &fov, float &magn, vector<button>& buttons,
bool &showMenu, WINDOW* &win, int &currBtn)
{
        static int prevMouseX = -1, prevMouseY = -1;
        MEVENT event;
        ch = getch();
        // Обработка мыши
    if (ch == KEY_MOUSE && getmouse(&event) == OK) {
        // Смещение мыши
        if (prevMouseX != -1 && prevMouseY != -1) {
            int deltaX = event.x - prevMouseX;
            int deltaY = event.y - prevMouseY;

            // Поворот камеры по горизонтали (мышь X)
            float angleY = deltaX * 0.01f;
            float oldDirX = cam.direction.x;
            cam.direction.x = cam.direction.x * cos(-angleY) - cam.direction.z * sin(-angleY);
            cam.direction.z = oldDirX * sin(-angleY) + cam.direction.z * cos(-angleY);

            // Наклон камеры по вертикали (мышь Y)
            float angleX = deltaY * 0.01f;
            cam.direction.y -= angleX;
        }
        prevMouseX = event.x;
        prevMouseY = event.y;
        // Обработка колесика мыши
        if (event.bstate & BUTTON5_PRESSED) {  // Прокрутка вверх
            if (magn > 1) magn--;
            } else if (event.bstate & BUTTON4_PRESSED) {  // Прокрутка вниз
            magn++;
            }
    }
    else {
        // Остальное
        if (ch == '+')
            w.globalScale += 0.1f * magn;
        else if (ch == '-')
            w.globalScale -= 0.1f * magn;
        else if (ch == 'w')
            cam.physics.force = vecAdd(cam.physics.force, multVecScal(normalize(cam.direction), magn)); //cam.position.z += 0.1f * magn;
        else if (ch == 's')
            cam.physics.force = vecSubtract(cam.physics.force, multVecScal(normalize(cam.direction), magn));
        else if (ch == 'z')
            cam.position.y -= 0.1f * magn;
        else if (ch == ' ')
            cam.position.y += 0.1f * magn;
        else if (ch == 'a')
            cam.position.x -= 0.1f * magn;
        else if (ch == 'd')
            cam.position.x += 0.1f * magn;
        else if (ch == ']')
            magn += 1.0f;
        else if (ch == '[')
            magn -= 1.0f;
        else if (ch == 'm') // Главное меню
        {
        showMenu = !showMenu;
        if (showMenu) {
                win = createWindow(buttons, 25, 50, -1, -1,
                                  COLOR_BLUE, COLOR_BLACK,
                                  "Main Menu", true);
                // Первоначальная отрисовка
                for (size_t i = 0; i < buttons.size(); ++i) {
                    updateButton(win, buttons[i], i == currBtn);
                }
            } else {
                //showMenu = false;
                delwin(win);
                win = nullptr;
                clear();
                refresh();
            }
        }

        if (showMenu && win) {
            ch = wgetch(win);

            switch (ch) {
                case KEY_UP:
                    if (currBtn > 0) {
                        currBtn--;
                    } else currBtn = buttons.size() - 1;
                    break;
                case KEY_DOWN:
                    if (currBtn < buttons.size() - 1) {
                        currBtn++;
                    } else currBtn = 0;
                    break;
                case KEY_LEFT:
                    if(buttons[currBtn].type == BUTTON) {
                        *buttons[currBtn].targetVar -= buttons[currBtn].step;
                    }
                    break;
                case KEY_RIGHT:
                    if(buttons[currBtn].type == BUTTON) {
                        *buttons[currBtn].targetVar += buttons[currBtn].step;
                    }
                    break;
                case KEY_ENTER: // Ввод для переключения чекбокса
                    if(buttons[currBtn].type == CHECKBOX) {
                        *buttons[currBtn].targetBool = !*buttons[currBtn].targetBool;
                    }
                case '\n': // Ввод для переключения чекбокса
                    if(buttons[currBtn].type == CHECKBOX) {
                        *buttons[currBtn].targetBool = !*buttons[currBtn].targetBool;
                    }
                    break;
                case 'q':
                    showMenu = false;
                    delwin(win);
                    win = nullptr;
                    clear();
                    refresh();
                    break;
            }

            // Проверка границ значений (только для BUTTON)
            if(buttons[currBtn].type == BUTTON) {
            button& btn = buttons[currBtn];
            if (*btn.targetVar < btn.minVal) *btn.targetVar = btn.minVal;
            if (*btn.targetVar > btn.maxVal) *btn.targetVar = btn.maxVal;
            }

            // Обновление отображения
            for (size_t i = 0; i < buttons.size(); ++i) {
                updateButton(win, buttons[i], i == currBtn);
            }
        }
            // Защита от слишком маленького масштаба
        if (w.globalScale < 0.01f)
            w.globalScale = 0.01f;
	    if(magn < 0)
	        magn = 0.0f;
    }
}

// Счет треугольников
size_t cntTris(const world& w)
{
    size_t t = 0;
    for (const auto& obj : w.objects)
    {
        t+= obj.model.triangles.size();
    }
    return t;
}

// Перевод в сферические
vec3d toSpherical(const float& elev, const float& azim, const float& magn) {
    // elev - угол подьема (0 - 90), azim - азимут (0 - 360), magn - радиус (магнитуда)
    float radE = elev * (M_PI / 180);
    float radA = azim * (M_PI / 180);
    return {
        magn * cos(radE) * cos(radA), // X = r * cos(radE) * cos(radA)
        magn * sin(radE),     // Y = r * sin(radE)
        magn * cos(radE) * sin(radA)  // Z = r * cos(radE) * sin(radA)
    };
}

// Красивый вывод выбора поворота
void arrow(const mesh& model, const float& elev, const float& azim, const float& posX, const float& posY, const float& posZ, const matrx4d &projMatrix)
{
    matrx4d rotMatrix = matMult(createRotationX(azim), createRotationY(elev));
    matrx4d viewMatrix = createViewMatrix(cam);
    matrx4d transMatrix = createTranslationMatrix(posX, posY, posZ);
    drawMesh(model, projMatrix, viewMatrix, 1.0f, rotMatrix, transMatrix, cam);
}

void updateCam(camera& cam)
{
    cam.position = cam.physics.pos;
}

// Начальный экран
void startScrn(world* w)
{
object obj;
bool start = false;
bool nObj = false;
bool put = false;
int currBtn = 0;
int ch;
float magn;
float elev;
float azim;
vector<button> buttons  = {
{CHECKBOX, static_cast<int>(LINES /  2) - 8, 10, .targetBool = &nObj, "Create New Object"}};
button cw = {CHECKBOX, static_cast<int>(LINES / 2) - 8, 12, .targetBool = &start, "Create New World"};

while(!start) {
      getmaxyx(stdscr, height, width);
	  if(w->objects.size() > 1  && !put)
	  {
	    buttons.push_back(cw);
	    put = true;
	  }
	  WINDOW* win = createWindow(buttons, height, width / 2, -1, 0,
                                  COLOR_WHITE, COLOR_BLACK,
                                  "Welcome to TraceIt!"
                                  , true);
                // Первоначальная отрисовка
                for (size_t i = 0; i < buttons.size(); ++i) {
                    updateButton(win, buttons[i], i == currBtn);
                }
                mvwprintw(win, 10, static_cast<int>(LINES /  2) + 13, " (In Total:  %i)", static_cast<int>(w->objects.size()));

        if (win) {
            ch = wgetch(win);

            switch (ch) {
                case KEY_UP:
                    if (currBtn > 0) {
                        currBtn--;
                    } else currBtn = buttons.size() - 1;
                    break;
                case KEY_DOWN:
                    if (currBtn < buttons.size() - 1) {
                        currBtn++;
                    } else currBtn = 0;
                    break;
                case KEY_LEFT:
                    if(buttons[currBtn].type == BUTTON) {
                        *buttons[currBtn].targetVar -= buttons[currBtn].step;
                    }
                    break;
                case KEY_RIGHT:
                    if(buttons[currBtn].type == BUTTON) {
                        *buttons[currBtn].targetVar += buttons[currBtn].step;
                    }
                    break;
                case KEY_ENTER: // Ввод для переключения чекбокса
                    if(buttons[currBtn].type == CHECKBOX) {
                        *buttons[currBtn].targetBool = !*buttons[currBtn].targetBool;
                    }
                    break;
                case '\n': // Ввод для переключения чекбокса
                    if(buttons[currBtn].type == CHECKBOX) {
                        *buttons[currBtn].targetBool = !*buttons[currBtn].targetBool;
                    }
                    break;
            }
	if(nObj)  {
	clear();
	refresh();
	printw("\nEnter Mass (kg): \n");
	try {
       		obj.mass = stof(consoleIn());
   		} catch (const std::invalid_argument& e) {
       		// Обработка ошибки
       		obj.mass = 1.0f;
     		printw("Not a number, defaulting to: %.2f", obj.mass);
  	    }
	printw("\nEnter Cd: \n");
	try {
                obj.Cd = stof(consoleIn());
        } catch (const std::invalid_argument& e) {
                // Обработка ошибки
                obj.Cd = 0.5f;
                printw("Not a number, defaulting to: %.2f", obj.Cd);
         }
	printw("\nVelocity in spherical, elevation (0 - 90): \n");
	try {
                elev = stof(consoleIn());
        } catch (const std::invalid_argument& e) {
                // Обработка ошибки
                elev = 0.0f;
                printw("Not a number, defaulting to: %.2f", elev);
         }
	printw("\nVelocity in spherical, azimuth (0 - 360): \n");
	try {
                azim = stof(consoleIn());
        } catch (const std::invalid_argument& e) {
                // Обработка ошибки
                azim = 0.0f;
                printw("Not a number, defaulting to: %.2f", azim);
           }
    printw("\nEnter velocity magnitude (m/s): \n");
	try {
                magn = stof(consoleIn());
        } catch (const std::invalid_argument& e) {
                // Обработка ошибки
                magn = 0.0f;
                printw("Not a number, defaulting to: %.2f", magn);
           }
    arrow(loadFromObjectFile("arrow.obj"), elev, azim, width - width / 4, height / 2, 0.0f, createProjMatrix(width / height, 90.0f));
    obj.vel = toSpherical(elev, azim, magn);
    printw("\nEnter Position (X): \n");
        try {
                obj.pos.x = stof(consoleIn());
            } catch (const std::invalid_argument& e) {
                // Обработка ошибки
                obj.pos.x = 0.0f;
                printw("Not a number, defaulting to: %.2f", obj.pos.x);
           }
    printw("\nEnter Position (Y): \n");
        try {
                obj.pos.y = stof(consoleIn());
            } catch (const std::invalid_argument& e) {
                // Обработка ошибки
                obj.pos.y = 0.0f;
                printw("Not a number, defaulting to: %.2f", obj.pos.y);
           }
    printw("\nEnter Position (Z): \n");
        try {
                obj.pos.z = stof(consoleIn());
                } catch (const std::invalid_argument& e) {
                // Обработка ошибки
                obj.pos.z = 0.0f;
                printw("Not a number, defaulting to: %.2f", obj.pos.z);
           }
    obj.transform = createTranslationMatrix(obj.pos.x, obj.pos.y, obj.pos.z);
 	obj.force = {0, 0, 0};
 	obj.area = 2.0f;
	obj.model = loadFromObjectFile("sphere.obj");
	// Инициализация коллайдера
    obj.col.center = obj.pos;
    obj.col.size = {2, 2, 2};
    obj.col.isSphere = false;
    obj.col.isStatic = false;
    obj.col.rest = 0.5;
    obj.stationary = false;
    w->objects.push_back(obj);
    nObj = false;
	clear();
	refresh();
	}
}

            // Проверка границ значений (только для BUTTON)
            if(buttons[currBtn].type == BUTTON) {
            button& btn = buttons[currBtn];
            if (*btn.targetVar < btn.minVal) *btn.targetVar = btn.minVal;
            if (*btn.targetVar > btn.maxVal) *btn.targetVar = btn.maxVal;
            }

            // Обновление отображения
            for (size_t i = 0; i < buttons.size(); ++i) {
                updateButton(win, buttons[i], i == currBtn);
            }
	}
	float ad;
	float g;
	clear();
	refresh();
	printw("\nEnter atmosphere density: \n");
	try {
		if((ad = stof(consoleIn())) <  0.0f) {
        	w->atmDensity = 1.125f;
        	}
        	else
       		w->atmDensity = ad;
                } catch (const std::invalid_argument& e) {
                w->atmDensity = 1.125;
                printw("Not a number, defaulting to: %.2f", w->atmDensity);
	}
	printw("\nEnter gravitational acceleration: \n");
	 try {
		if((g = stof(consoleIn())) <  0.0f) {
        	w->g = -g;
        	clear();
        	refresh();
        	}
        	else
       		w->atmDensity = g;
                } catch (const std::invalid_argument& e) {
                printw("Not a number, defaulting to: %.2f", w->g);
                clear();
	            refresh();
	}
}

// Главная функция
int main()
{
    // Инициализация ncurses
    initscr();
    cbreak();               // Без enter
    noecho();               // Отключение вывода
    keypad(stdscr, TRUE);   // Включение клавиш по типу стрелок
    nodelay(stdscr, TRUE);  // Включаем неблокирующий ввод
    curs_set(0);             // Скрываем курсор
    start_color();          // Включение цвета (пока только для меню)
    // Мышь
    mousemask(ALL_MOUSE_EVENTS | REPORT_MOUSE_POSITION, NULL);
    mouseinterval(0);  // Отключаем задержку для мыши
    scrollok(stdscr, TRUE);  // Включаем скроллинг
    printf("\033[?1003h\n"); // Для iterm?

    // Важные переменные
    int frames = 0;
    int currBtn = 0;
    float rotAngle = 0.0f;
    float magn = 1.0f;
    float fov = 90.0;
    float objIndex = 0;
    bool isRunning = true;
    bool debug = true;
    bool enablePhysics = false;
    bool enableTriCnt = false;
    bool showMenu = false;
    bool disableFilling = true;
    bool objView = false;
    bool countTillFallen = false; // Доделать
    cam.physics = {{0.0f, 0.0f, -20.0f},
    {0.0f, 0.0f, 20.0f},
    {0.0f, 0.0f, 0.0f},
    75.0f, createRectMesh(0.4f, 1.75f, 0.25), createTranslationMatrix(0, 0, -20),
    0.7f, 0.1f, false};
    //cam.physics.trajWritten = true;
    cam.physics.col = {cam.position, {0.4f, 1.75f, 0.25}, 0.0f, false, false, 0.1};

    matrx4d floPos = createTranslationMatrix(0, -30, 0);

    // Создание мира
    mesh m = createRectMesh(2000, 5, 2000);
    world w;
    w.globalScale = 1.0f;
    object floor = {{0, -30, 0}, {0, 0, 0}, {0, 0, 0}, 1000.0f, m, floPos, 10.0f, 4.0f, true};
    floor.col.center = floor.pos;
    floor.col.size = {2000, 5, 2000};
    floor.col.isSphere = false;
    floor.col.isStatic = true;
    floor.col.rest = 1.0f;
    w.objects.push_back(floor);
    w.objects.push_back(cam.physics);
    startScrn(&w);

    // Создаем таблицу координат
    ofstream coords;

    // Матрицы
    matrx4d scaleMatrix;
    matrx4d rotMatrix;
    matrx4d rotX;
    matrx4d rotY;
    matrx4d rotZ;

    // Кнопки главного меню
    vector<button> mainM = {
        // Кнопки
        {BUTTON, 10, 3, .targetVar = &fov, "FOV", 50, 170, 5},
        {BUTTON, 10, 5, .targetVar = &w.atmDensity, "Atm. Density", 0, 50, 0.1},
        {BUTTON, 10, 7, .targetVar = &objIndex, "Object", 0, static_cast<int>(w.objects.size()), 1},
        {BUTTON, 10, 9, .targetVar = &magn, "Magnitude", 0, 500, 1},
        {BUTTON, 10, 11, .targetVar = &w.objects[objIndex].Cd, "Obj. Cd", 0, 50, 1},

        // Флажки
        {CHECKBOX, 10, 13, .targetBool = &enablePhysics, "Enable Physics"},
        {CHECKBOX, 10, 15, .targetBool = &enableTriCnt, "Count All Triangles"},
        {CHECKBOX, 10, 17, .targetBool = &debug, "Debug Mode"},
        {BUTTON, 10, 19, .targetVar = &w.objects[objIndex].vel.x, "Vel X", 0, 50, 1},
        {CHECKBOX, 10, 21, .targetBool = &objView, "Object PoW"},
        {CHECKBOX, 10, 24, .targetBool = &isRunning, "[QUIT]"}
        };

    WINDOW* win = nullptr;

    // Главный цикл
    while (isRunning)
    {
        // Получение и вывод размера консоли
        getmaxyx(stdscr, height, width);
        float aspect = static_cast<float>(height) / static_cast<float>(width);
        if(!showMenu) {
            clear(); // Очистка экрана
            rotAngle += 0.00f;
            rotX = createRotationX(rotAngle);
            rotY = createRotationY(rotAngle);
	        rotZ = createRotationZ(rotAngle);

            rotMatrix = matMult(rotX, rotY);
	        //rotMatrix = matMult(rotMatrix, rotZ); // Комбинирование поворотов
            // Создание матрицы проекции
            w.projMatrix = createProjMatrix(aspect, fov);
            // Создание матрицы вида для камеры
	        w.viewMatrix = createViewMatrix(cam);
            // Создание матрицы масштаба
	        scaleMatrix = createScaleMatrix(w.globalScale, w.globalScale, w.globalScale);

	        if(enablePhysics) updatePhysics(w, coords);
	        drawWorld(w, rotMatrix, cam);
	        updateCam(cam);
	        if(enableTriCnt) mvprintw(height - 3, 0, "All Triangles: %.i", cntTris(w));

            if(objView) {
                cam.direction = normalize(w.objects[objIndex].vel);
                cam.position = w.objects[objIndex].pos;
            }

            if(debug){
                move(4, 0);
                printw("Scale: %f\n", w.globalScale);
                move(5, 0);
                printw("Camera position: (x%.2f, y%.2f, z%.2f tiltY %.2f, tiltX %.2f) \n", cam.position.x, cam.position.y, cam.position.z, cam.direction.y, cam.direction.x);
	            move(6, 0);
	            printw("Frame #%i, Time: %.2f", frames, (frames * 0.0166));
	            move(7, 0);
	            printw("Obj. Pos. x %.2f, y %.2f, z %.2f", w.objects[objIndex].pos.x, w.objects[objIndex].pos.y, w.objects[objIndex].pos.z);
	            move(8, 0);
	            printw("Obj. Vel. x %.2f, y %.2f, z %.2f", w.objects[objIndex].vel.x, w.objects[objIndex].vel.y, w.objects[objIndex].vel.z);
	            move(9, 0);
	            printw("Terminal Dimensions: Columns %i, Rows %i,", width, height);
        }
    }
        // Обработка ввода
        int ch;
        handleInput(ch, cam, w, fov, magn, mainM, showMenu, win, currBtn);

        refresh(); // Обновление экрана

        // Примерно 60 к/c (без лага)
        frames++;
        usleep(16600);
    }

    // Завершение работы ncurses
    endwin();
    return 0;
}



