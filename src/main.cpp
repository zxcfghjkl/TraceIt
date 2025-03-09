// Для ?
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

// Треугольник
struct triangle
{
    vec3d p[3];
};

// Сетка
struct mesh {
    vector<triangle> triangles; // Массив треугольников, составляющих сетку
};


// Камера
struct camera
{
    vec3d position;  // Позиция камеры
    vec3d direction; // Направление взгляда
    vec3d up;        // Вектор "вверх" камеры
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
    vector<vec3d> trajectory; // История позиций
    int maxTrajectoryPoints = 300; // Макс. точек траектории
    bool trajWritten = false;
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
camera cam = {
    {0.0f, 0.0f, -140.0f},  // Позиция камеры (сдвинута назад)
    {0.0f, 0.0f, 1.0f},   // Направление взгляда (смотрит по оси Z)
    {0.0f, -1.0f, 0.0f}    // Вектор "вверх" (по оси Y)
};

// Функция, использующая алгоритм Брезенхема для рисования линий
void drawLine(int x1, int y1, int x2, int y2)
{
    int dx = abs(x2 - x1);       // Разница по X
    int dy = abs(y2 - y1);       // Разница по Y
    int sx = (x1 < x2) ? 1  : -1; // Направление движения по X
    int sy = (y1 < y2) ? 1 : -1; // Направление движения по Y
    int err = dx - dy;           // Начальная ошибка

    while (true)
    {
        mvaddch(y1, x1, '*'); // Рисуем символ '*' в текущей позиции
        if (x1 == x2 && y1 == y2)
            break; // Если достигли конца линии, выходим из цикла

        int e2 = 2 * err;
        if (e2 > -dy)
        { // Корректируем ошибку для X
            err -= dy;
            x1 += sx;
        }
        if (e2 < dx)
        { // Корректируем ошибку для Y
            err += dx;
            y1 += sy;
        }
    }
}

// Точно не скопировано у OneLoneCoder
mesh loadFromObjectFile(string sFilename)
	{
	    mesh outputMesh;

		ifstream f(sFilename);

		// Локальный кеш вершин
		vector<vec3d> verts;

		while (!f.eof())
		{
			char line[128];
			f.getline(line, 128);

			strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				vec3d v;
				s >> junk >> v.x >> v.y >> v.z;
				verts.push_back(v);
			}

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
// Нужно для определения видимости
vec3d calculateNormal(const triangle &tri)
{
    vec3d edge1 = vecSubtract(tri.p[1], tri.p[0]);
    vec3d edge2 = vecSubtract(tri.p[2], tri.p[0]);
    return normalize(cross(edge1, edge2));
}

// А вот и она
bool isTriangleVisible(const triangle &tri, const camera& cam)
{
    vec3d normal = calculateNormal(tri);
    vec3d toTriangle = vecSubtract(tri.p[0], cam.direction);

    return dot(normal, toTriangle) > 0;// Если скалярное произведение отрицательное, треугольник направлен в сторону камеры
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
    for (int i = 0; i < 3; i++) {
        // Правильный порядок: Масштаб > Поворот > Перемещение > Вид > Проекция
        vec3d point = tri.p[i];
        point = matVecMult(scaleMatrix, point);
        point = matVecMult(rotMatrix, point);
    	point = matVecMult(transMatrix, point);
        point = matVecMult(viewMatrix, point);
        point = matVecMult(projMatrix, point);
        tri.p[i] = scaleToScreen(point);
    }

    drawLine(tri.p[0].x, tri.p[0].y, tri.p[1].x, tri.p[1].y);
    drawLine(tri.p[1].x, tri.p[1].y, tri.p[2].x, tri.p[2].y);
    drawLine(tri.p[2].x, tri.p[2].y, tri.p[0].x, tri.p[0].y);

    move(10, 0);

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
        matrx4d projMatrix = {{ // 2.333 - соотношение сторон символа
          {fovRad * aspect * 2.333f, 0.0f, 0.0f, 0.0f},
          {0.0f,            fovRad,    0.0f,                 0.0f},
          {0.0f,            0.0f,      (far + near)/(far - near), 1.0f},
          {0.0f,            0.0f,      (-2.0f * far * near)/(far - near), 0.0f}
          }};

        return projMatrix;
}

// Рисуем сетку
void drawMesh(const mesh& inputMesh, const matrx4d& projMatrix, const matrx4d& viewMatrix, float globalScale,
 		const matrx4d& rotMatrix, const matrx4d& transMatrix, const camera& cam) {
    // Создаём матрицу масштаба
    matrx4d scaleMatrix = createScaleMatrix(globalScale, globalScale, globalScale);

    // Проходимся по всем треугольникам в сетке
    for (const triangle& trian : inputMesh.triangles) {

        // Проверяем видимость треугольника
        if (isTriangleVisible(trian, cam)) {
            drawTriangle(trian, projMatrix, viewMatrix, scaleMatrix, rotMatrix, transMatrix); // Отрисовка видимого треугольника
        }
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

object createCubeObject(const float &scale, const vec3d &pos,
                        const vec3d &vel, const vec3d &force,
                        const float &mass, const float &area, const float &Cd) {
    mesh cubeMesh;

    // Масштабирование вершины
    auto scaleVertex = [scale](vec3d v) -> vec3d {
        return {v.x * scale, v.y * scale, v.z * scale};
    };

    // Вершины куба
    vec3d vertices[8] = {
        {-0.5f, -0.5f, -0.5f}, // 0
        {0.5f, -0.5f, -0.5f},  // 1
        {0.5f, 0.5f, -0.5f},   // 2
        {-0.5f, 0.5f, -0.5f},  // 3
        {-0.5f, -0.5f, 0.5f},  // 4
        {0.5f, -0.5f, 0.5f},   // 5
        {0.5f, 0.5f, 0.5f},    // 6
        {-0.5f, 0.5f, 0.5f}    // 7
    };

    // Определяем треугольники для каждой грани куба
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

    // Применяем масштаб
    for(auto& tri : cubeTriangles) {
        for(auto& v : tri.p) {
            v.x *= scale;
            v.y *= scale;
            v.z *= scale;
        }
    }

    // Добавляем треугольники в mesh
    for (const auto& tri : cubeTriangles) { // Для каждого треугольника в эррее "cubeTriangles"
        cubeMesh.triangles.push_back(tri);
    }

    // Создаем обьект
    object cubeObj;
    cubeObj.pos = pos;
    cubeObj.transform = createTranslationMatrix(pos.x, pos.y, pos.z);
    cubeObj.vel = vel;
    cubeObj.force = force;
    cubeObj.mass = mass;
    cubeObj.Cd = Cd;
    cubeObj.area = area;
    cubeObj.model.triangles.assign(begin(cubeTriangles), end(cubeTriangles));

    return cubeObj;
}

// Обновление матрицы трансляции в обьекте
void updateObjectTransform(object& obj) {
    obj.transform = createTranslationMatrix(obj.pos.x, obj.pos.y, obj.pos.z);
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
        ch = getch();
        if (ch == '+')       // Если нажата клавиша '+'
            w.globalScale += 0.1f * magn; // Увеличить масштаб
        else if (ch == '-')  // Если нажата клавиша '-'
            w.globalScale -= 0.1f * magn; // Уменьшить масштаб
        else if (ch == KEY_UP) // Перемещение камеры вперед
            cam.position.z += 0.1f * magn;
        else if (ch == KEY_DOWN) // Перемещение камеры назад
            cam.position.z -= 0.1f * magn;
        else if (ch == 'd') // Перемещение камеры вверх
            cam.position.y -= 0.1f * magn;
        else if (ch == 'a') // Перемещение камеры вниз
            cam.position.y += 0.1f * magn;
        else if (ch == 's') // Наклон вверх
            cam.direction.y -= 0.1f;
        else if (ch == 'w') // Наклон вниз
            cam.direction.y += 0.1f;
        else if (ch == 'z') // Влево
            cam.position.x -= 0.1f * magn;
        else if (ch == 'x') // Вправо
            cam.position.x += 0.1f * magn;
        else if (ch == KEY_LEFT) // Поворот камеры влево
        {
            float angle = 0.1f;
            cam.direction.x = cam.direction.x * cos(angle) - cam.direction.z * sin(angle);
            cam.direction.z = cam.direction.x * sin(angle) + cam.direction.z * cos(angle);
        }
        else if (ch == KEY_RIGHT) // Поворот камеры вправо
        {
            float angle = -0.1f;
            cam.direction.x = cam.direction.x * cos(angle) - cam.direction.z * sin(angle);
            cam.direction.z = cam.direction.x * sin(angle) + cam.direction.z * cos(angle);
        }else if (ch == ']') // Увеличение магнитуды
            magn += 1.0f;
	    else if (ch == '[') // Уменьшение магнитуды
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

// Начальный экран
void startScrn(world* w)
{
object obj;
bool start = false;
bool nObj = false;
int currBtn = 0;
int ch;
bool put = false;
vector<button> buttons  = {
{CHECKBOX, static_cast<int>((LINES - 5) /  2), 13, .targetBool = &nObj, "Create New Object"}};
button cw = {CHECKBOX, static_cast<int>((LINES - 2) / 2), 15, .targetBool = &start, "Create New World"};

while(!start) {
          getmaxyx(stdscr, height, width);
	  if(w->objects.size() > 0  && !put)
	  {
	  buttons.push_back(cw);
	  put = true;
	  }
	  WINDOW* win = createWindow(buttons, height - 1, width - 4, -1, -1,
                                  COLOR_WHITE, COLOR_BLACK,
                                  " _____                     _____ _     _\n     /__   \\_ __ __ _  ___ ___  \\_   \\ |_  / \\\n/ /\\/ '__/ _` |/ __/ _ \\  / /\\/ __|/  /\n    / /  | | | (_| | (_|  __/\\/ /_ | |_/\\_/\n\\/   |_|  \\__,_|\\___\\___\\____/  \\__\\/"
                                  , true);
                // Первоначальная отрисовка
                for (size_t i = 0; i < buttons.size(); ++i) {
                    updateButton(win, buttons[i], i == currBtn);
                }
                wprintw(win, " (Total objects:  %i)", static_cast<int>(w->objects.size()));

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
	printw("\nEnter Mass: \n");
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
	printw("\nEnter Velocity (X): \n");
	try {
                obj.vel.x = stof(consoleIn());
                } catch (const std::invalid_argument& e) {
                // Обработка ошибки
                obj.vel.x = 0.0f;
                printw("Not a number, defaulting to: %.2f", obj.vel.x);
           }
	printw("\nEnter Velocity (Y): \n");
	try {
                obj.vel.y = stof(consoleIn());
                } catch (const std::invalid_argument& e) {
                // Обработка ошибки
                obj.vel.y = 0.0f;
                printw("Not a number, defaulting to: %.2f", obj.vel.y);
           }
	printw("\nEnter Velocity (Z): \n");
        try {
                obj.vel.z = stof(consoleIn());
                } catch (const std::invalid_argument& e) {
                // Обработка ошибки
                obj.vel.z = 0.0f;
                printw("Not a number, defaulting to: %.2f", obj.vel.z);
           }
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
 	obj.force = {0, 0, 0};
	obj.model = loadFromObjectFile("teapot.obj");
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
        	clear();
        	refresh();
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
	}
}

// Главная функция
int main()
{
    // Инициализация ncurses
    initscr();
    cbreak();
    noecho();               // Отключение вbвода
    keypad(stdscr, TRUE);   // Включение клавиш по типу стрелок
    nodelay(stdscr, TRUE);  // Включаем неблокирующий ввод
    curs_set(0);             // Скрываем курсор
    start_color();          // Включение цвета (пока только для меню)

    // Создание мира
    world w;
    startScrn(&w);

    // Создаем таблицу координат
    ofstream coords;
    bool enablePhysics = false;
    bool enableTriCnt = true;

    //for(int i = 0; i < 2; i++){
    //vec3d objPos = {static_cast<float>(-i * 6), 0, 0};
    //vec3d objVel = {0, 0, 0};
    //vec3d objForce = {0, 0, 0};
    //object cube = createCubeObject(0.5f, objPos, objVel, objForce, i * 1.0f, i * 1.09f, i * 0.1f);
    //cube.model = loadFromObjectFile("teapot.obj");
    //w.objects.push_back(cube);
    //}

    // Коеффициент масштаба и плотности атмосферы
    w.globalScale = 1.0f;
    //w.atmDensity = 1.225;


    float rotAngle = 0.0f;
    float magn = 1.0f;
    int frames = 0;
    float fov = 90.0;
    float objIndex = 0;
    bool isRunning = true;
    bool debug = false;

    // Кнопки главного меню
    vector<button> mainM = {
        // Кнопки
        {BUTTON, 10, 3, .targetVar = &fov, "FOV", 60, 120, 5},
        {BUTTON, 10, 5, .targetVar = &w.atmDensity, "Atm. Density", 0, 50, 0.1},
        {BUTTON, 10, 7, .targetVar = &objIndex, "Object", 0, static_cast<int>(w.objects.size()), 1},
        {BUTTON, 10, 9, .targetVar = &magn, "Magnitude", 0, 500, 1},
        {BUTTON, 10, 11, .targetVar = &w.objects[objIndex].Cd, "Obj. Cd", 0, 50, 1},

        // Флажки
        {CHECKBOX, 10, 13, .targetBool = &enablePhysics, "Enable Physics"},
        {CHECKBOX, 10, 15, .targetBool = &enableTriCnt, "Count All Triangles"},
        {CHECKBOX, 10, 17, .targetBool = &debug, "Debug Mode"},
        {BUTTON, 10, 19, .targetVar = &w.objects[objIndex].vel.x, "Vel X", 0, 50, 1},
        {CHECKBOX, 10, 24, .targetBool = &isRunning, "[QUIT]"}
        };

    WINDOW* win = nullptr;
    bool showMenu = false;
    int currBtn = 0;
    // Главный цикл
    while (isRunning)
    {
        if(!showMenu) {
        // Получение и вывод размера консоли
        getmaxyx(stdscr, height, width);
        float aspect = static_cast<float>(height) / static_cast<float>(width);

        rotAngle += 0.02f;
        matrx4d rotX = createRotationX(rotAngle);
        matrx4d rotY = createRotationY(rotAngle);
	      matrx4d rotZ = createRotationZ(rotAngle);

        matrx4d rotMatrix = matMult(rotX, rotY);
	    rotMatrix = matMult(rotMatrix, rotZ); // Комбинирование поворотов

        clear(); // Очистка экрана

        // Создание матрицы проекции
        w.projMatrix = createProjMatrix(aspect, fov);

        // Создание матрицы вида для камеры
	    w.viewMatrix = createViewMatrix(cam);

	      // Создание матрицы масштаба
	    matrx4d scaleMatrix = createScaleMatrix(w.globalScale, w.globalScale, w.globalScale);

	    if(enablePhysics) updatePhysics(w, coords);
	    drawWorld(w, rotMatrix, cam);
	    if(enableTriCnt) mvprintw(height - 3, 0, "All Triangles: %.i", cntTris(w));

        if(debug){
            move(4, 0);
            printw("Scale: %f\n", w.globalScale);
            move(5, 0);
            printw("Camera position: (x%.2f, y%.2f, z%.2f tilt%.2f) \n", cam.position.x, cam.position.y, cam.position.z, cam.direction.y);

	        move(6, 0);
	        printw("Frame #%i", frames);

	        move(7, 0);
	        printw("Obj. Pos. x %.2f, y %.2f, z %.2f", w.objects[objIndex].pos.x, w.objects[objIndex].pos.y, w.objects[objIndex].pos.z);

	        move(8, 0);
	        printw("Obj. Vel. x %.2f, y %.2f, z %.2f", w.objects[objIndex].vel.x, w.objects[objIndex].vel.y, w.objects[objIndex].vel.z);
        }
    }
            // Обработка ввода
        int ch;
        handleInput(ch, cam, w, fov, magn, mainM, showMenu, win, currBtn);

        refresh(); // Обновление экрана

            // Примерно 60 к/c
            frames++;
            usleep(16666);
    }

        // Завершение работы ncurses
    endwin();
    return 0;
}



