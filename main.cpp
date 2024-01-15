#include <iostream>
#include "rayos.h"
#include <fstream>
#include <sstream>

using namespace std;

//=======================================
//VariablesGlobales
bool salaIniciada; //si la habitación ha sido cargada .
room r;          // Instancia para representar la habitación.
MatrizResultados mTE; //matriz que guarda el ID del triángulo en el que rebota el rayo
MatrizResultados mE; //matriz de energia residual.
int NumTri = 0; //Almacena el número total de triángulos.
double N_RAYOS = 12; //Almacena el numero de rayos
int energiaFuente = 100; //Energia de la fuente
source s;           // Instancia de la fuente.
float alfa = 0.2; //Defincion del coeficiente alfa
float delta = 0.2; //Defincion del coeficiente alfa
reflexion* reflexiones = NULL; //Reflecciones de cada rayo
point o; //Punto de origen de la Fuente


void crearSala();
void calcularEnergia(double x, double y, double z);
void guardarMatrizCSV(MatrizResultados matriz, const std::string& nombreArchivo);

int main()
{
    //Se carga la sala
    crearSala();

    // Declaración de variables
    double val_X, val_Y, val_Z;

    // Mensaje para el usuario
    std::cout << "Ingrese la posicion x: ";

    // Lectura del primer valor desde la consola
    std::cin >> val_X;

    // Mensaje para el usuario
    std::cout << "Ingrese la posicion y: ";

    // Lectura del segundo valor desde la consola
    std::cin >> val_Y;

    // Mensaje para el usuario
    std::cout << "Ingrese la posicion z: ";

    // Lectura del tercer valor desde la consola
    std::cin >> val_Z;

    calcularEnergia(val_X, val_Y, val_Z);


    return 0;
}

//===================================Crear Sala=============================================================================================================
// Función para calculos al cargar la habitación:
void crearSala()
{
    if (!salaIniciada)     // Verifica si la habitación ya ha sido cargada anteriormente.
    {
        r.NewPlanes(6);// Genearra 6 planos
        //square back
        r.p[0].NewPoints(4);// Gnererar los 4 puntos
        r.p[0].p[0].x = -2.0f;
        r.p[0].p[0].y = 2.0f;
        r.p[0].p[0].z = 2.0f;
        r.p[0].p[1].x = -2.0f;
        r.p[0].p[1].y = -2.0f;
        r.p[0].p[1].z = 2.0f;
        r.p[0].p[2].x = -2.0f;
        r.p[0].p[2].y = -2.0f;
        r.p[0].p[2].z = -2.0f;
        r.p[0].p[3].x = -2.0f;
        r.p[0].p[3].y = 2.0f;
        r.p[0].p[3].z = -2.0f;
        r.p[0].PointGenTriangle();
        //square front
        r.p[1].NewPoints(4);// Gnererar los 4 puntos
        r.p[1].p[0].x = 2.0f;
        r.p[1].p[0].y = 2.0f;
        r.p[1].p[0].z = 2.0f;
        r.p[1].p[1].x = 2.0f;
        r.p[1].p[1].y = 2.0f;
        r.p[1].p[1].z = -2.0f;
        r.p[1].p[2].x = 2.0f;
        r.p[1].p[2].y = -2.0f;
        r.p[1].p[2].z = -2.0f;
        r.p[1].p[3].x = 2.0f;
        r.p[1].p[3].y = -2.0f;
        r.p[1].p[3].z = 2.0f;
        r.p[1].PointGenTriangle();
        //square left
        r.p[2].NewPoints(4);
        r.p[2].p[0].x = -2.0f;
        r.p[2].p[0].y = -2.0f;
        r.p[2].p[0].z = 2.0f;
        r.p[2].p[1].x = 2.0f;
        r.p[2].p[1].y = -2.0f;
        r.p[2].p[1].z = 2.0f;
        r.p[2].p[2].x = 2.0f;
        r.p[2].p[2].y = -2.0f;
        r.p[2].p[2].z = -2.0f;
        r.p[2].p[3].x = -2.0f;
        r.p[2].p[3].y = -2.0f;
        r.p[2].p[3].z = -2.0f;
        r.p[2].PointGenTriangle();
        //square right
        r.p[3].NewPoints(4);// Gnererar los 4 puntos
        r.p[3].p[0].x = 2.0f;
        r.p[3].p[0].y = 2.0f;
        r.p[3].p[0].z = 2.0f;
        r.p[3].p[1].x = -2.0f;
        r.p[3].p[1].y = 2.0f;
        r.p[3].p[1].z = 2.0f;
        r.p[3].p[2].x = -2.0f;
        r.p[3].p[2].y = 2.0f;
        r.p[3].p[2].z = -2.0f;
        r.p[3].p[3].x = 2.0f;
        r.p[3].p[3].y = 2.0f;
        r.p[3].p[3].z = -2.0f;
        r.p[3].PointGenTriangle();
        //square top
        r.p[4].NewPoints(4);
        r.p[4].p[0].x = -2.0f;
        r.p[4].p[0].y = -2.0f;
        r.p[4].p[0].z = 2.0f;
        r.p[4].p[1].x = -2.0f;
        r.p[4].p[1].y = 2.0f;
        r.p[4].p[1].z = 2.0f;
        r.p[4].p[2].x = 2.0f;
        r.p[4].p[2].y = 2.0f;
        r.p[4].p[2].z = 2.0f;
        r.p[4].p[3].x = 2.0f;
        r.p[4].p[3].y = -2.0f;
        r.p[4].p[3].z = 2.0f;
        r.p[4].PointGenTriangle();
        //square bottom
        r.p[5].NewPoints(4);
        r.p[5].p[0].x = -2.0f;
        r.p[5].p[0].y = 2.0f;
        r.p[5].p[0].z = -2.0f;
        r.p[5].p[1].x = -2.0f;
        r.p[5].p[1].y = -2.0f;
        r.p[5].p[1].z = -2.0f;
        r.p[5].p[2].x = 2.0f;
        r.p[5].p[2].y = -2.0f;
        r.p[5].p[2].z = -2.0f;
        r.p[5].p[3].x = 2.0f;
        r.p[5].p[3].y = 2.0f;
        r.p[5].p[3].z = -2.0f;
        r.p[5].PointGenTriangle();

        int cont_t = 0; //contador del numero de triangulos
        // Loop para calcular el centroide de todos los triángulos de la sala
        for (int i = 0; i < r.NP; i++)                   // Recorre los planos de la sala.
        {
            r.p[i].n = r.p[i].NormalPlano();       // Calcula la normal del plano.
            for (int j = 0; j < r.p[i].NT; j++)          // Recorre los triángulos del plano.
            {
                r.p[i].t[j].CalculateProjection();       // Calcula la proyección del triángulo.
                r.p[i].t[j].Centroid();                  // Calcula el centroide baricentro del triángulo.
                r.p[i].t[j].ID = cont_t;                 // Asigna un ID al triángulo.
                cont_t++;                                // Incrementa el contador de triángulos.
            }
        }
        NumTri = cont_t;  // Asigna el número total de triángulos.
        salaIniciada = true; // Indica que la sala ha sido cargada exitosamente.

    }
}

void calcularEnergia(double x, double y, double z)
{

    if (salaIniciada)
    {

        // Coordenadas de la fuente
        o.x = x;
        o.y = y;
        o.z = z;

        s.eF = energiaFuente; //Energia de la fuente

        s.createRays(N_RAYOS);  //Numero de rayos

        double eneRayo = s.eF / s.NRAYS; // Energía Inicial de cada rayo

        cout << "\n\nEnergia inicial de cada rayo: " << eneRayo << "\n";

        reflexiones = r.RayTracing(o, s.Rays, s.NRAYS);

        mE.Init(N_RAYOS, NUM_REBOTES); //inicializacion de la matriz de energia
        mTE.Init(N_RAYOS, NUM_REBOTES); //inicializacion de la matriz de ID de triangulos en el que rebota el rayo

        // Calculo Energia residual
        for (int rayo = 0; rayo < s.NRAYS; rayo++)
        {
            double eneResidual = eneRayo;
            for (int re = 0; re < reflexiones[rayo].N; re++)
            {
                int tri = reflexiones[rayo].idTriangle[re];
                eneResidual = eneResidual * (1 - alfa) * (1 - delta); //Energia residual en cada rebote
                mE.A[rayo][re] = eneResidual;
                mTE.A[rayo][re] = tri;

            }
        }

        std::cout << "\nMatriz de energia residual en cada rebote:" << std::endl;
        for (int i = 0; i < N_RAYOS; i++)
        {
            for (int j = 0; j < NUM_REBOTES; j++)
            {
                std::cout << mE.A[i][j] << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "\nMatriz ID de triangulos donde rebota el rayo:" << std::endl;
        for (int i = 0; i < N_RAYOS; i++)
        {
            for (int j = 0; j < NUM_REBOTES; j++)
            {
                std::cout << mTE.A[i][j] << " ";
            }
            std::cout << std::endl;
        }
//-------------------------------Guardar en Archivo CSV----------------------------------
        // Nombre del archivo CSV
        std::string nombreArchivoE = "matrizEnergiaResidual.csv";

        // Llamar a la función para guardar la matriz en el archivo CSV
        guardarMatrizCSV(mE, nombreArchivoE);

        std::cout << "\nMatriz guardada en " << nombreArchivoE << std::endl;

        // Nombre del archivo CSV
        std::string nombreArchivoT = "matrizIDtriangulosDondeRebotaelRayo.csv";

        // Llamar a la función para guardar la matriz en el archivo CSV
        guardarMatrizCSV(mTE, nombreArchivoT);

        std::cout << "Matriz guardada en " << nombreArchivoT << std::endl;



    }
}


void guardarMatrizCSV(MatrizResultados matriz, const std::string& nombreArchivo)
{
    std::ofstream archivo(nombreArchivo);

    if (!archivo.is_open())
    {
        std::cerr << "Error al abrir el archivo." << std::endl;
        return;
    }

    // Escribir los datos de la matriz en el archivo CSV
    for (int i = 0; i < N_RAYOS; i++)
    {
        for (int j = 0; j < NUM_REBOTES; j++)
        {
            archivo << matriz.A[i][j] << ";";
        }
        archivo << "\n";  // Nueva línea después de cada fila


    }
    archivo.close();
}
