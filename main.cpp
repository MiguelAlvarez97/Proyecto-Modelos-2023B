#include <iostream>
#include "rayos.h"
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

//=======================================
//Variables Globales
bool salaIniciada;      //si la habitación ha sido cargada .
room r;                 // Instancia para representar la habitación.
int tM = 1000;          //tiempo de simulacion en milisegundos
MatrizResultados mTE;   //matriz que guarda el ID del triángulo en el que rebota el rayo
MatrizResultados mE;    //matriz de energia residual.
int NumTri = 0;         //Almacena el número total de triángulos.
double N_RAYOS = 12;    //Almacena el numero de rayos
int energiaFuente = 120;//Energia de la fuente
source s;               // Instancia de la fuente.
float alfa = 0.2;       //Defincion del coeficiente alfa
float delta = 0.2;      //Defincion del coeficiente delta
reflexion* reflexiones = NULL; //Reflecciones de cada rayo
point o;                //Punto de origen de la Fuente
int cortes = 4;         //Número de cortes por cada plano
point** arrayRec;       //Array Recorrido
MatrizResultados mD;    //Matriz que almacene la distancia que existe entre los centros de los diferentes triángulos de la sala.
                        //Está matriz tiene una dimensión de 𝑛 × 𝑛 donde 𝑛 es el número de triángulos (NumTri) de la sala.
MatrizResultados mTV;   //matriz que almacena los tiempos de vuelo entre baricentros en milisegundos que demorarían las reflexiones difusas
                        //para llegar de un centroide a otro(entre todos los triángulos que apliquen, es decir que sean visibles)
                        // La matriz tendrá una dimensión de 𝑛 × n
MatrizResultados mAS; //matriz angulos solidos
MatrizResultados mPE; //matriz porcentage de energia. La matriz tendrá una dimensión de 𝑛 × n
MatrizResultados mET; //matriz porcentage de energia. (espacio/tiempo)
MatrizResultados mDR; //Matriz Distancia Receptor
MatrizResultados mTVR; //Matriz Tiempo de Vuelo  Receptor
MatrizResultados mASR; //matriz angulos solidos Receptor
MatrizResultados mPER; //matriz porcentage de energia Receptor
MatrizResultados mER; //matriz de energía residual

const int NumReceptores = 27; // Número de receptores

void crearSala();
void calcularEnergia(double x, double y, double z);
void guardarMatrizCSV(MatrizResultados matriz, const string& nombreArchivo);

int main()
{
    //Se carga la sala
    crearSala();

    // Declaración de variables
    double val_X, val_Y, val_Z;

    // Mensaje para el usuario
    cout << " \nIngrese la posición de la Fuente: " << endl;
    cout << "Ingrese la posiciones de la Fuente en x: ";

    // Lectura del primer valor desde la consola
    cin >> val_X;

    // Mensaje para el usuario
    cout << "Ingrese la posicion de la Fuente en y: ";

    // Lectura del segundo valor desde la consola
    cin >> val_Y;

    // Mensaje para el usuario
    cout << "Ingrese la posicion de la Fuente en z: ";

    // Lectura del tercer valor desde la consola
    cin >> val_Z;

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
        r.p[0].MoreTriangle(cortes);
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
        r.p[1].MoreTriangle(cortes);
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
        r.p[2].MoreTriangle(cortes);
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
        r.p[3].MoreTriangle(cortes);
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
        r.p[4].MoreTriangle(cortes);
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
        r.p[5].MoreTriangle(cortes);

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

        //salaIniciada = true; // Indica que la sala ha sido cargada exitosamente.
        // CREACIÓN DE RECEPTORES
        r.NewReceptor(NumReceptores);
        int cont_rec = 0;
        for (int i = -1; i < 2; i++) {
            for (int j = -1; j < 2; j++) {
                for (int k = -1; k < 2; k++) {
                    r.r[cont_rec].p.x = float(i);
                    r.r[cont_rec].p.y = float(j);
                    r.r[cont_rec].p.z = float(k);
                    cont_rec++;
                }
            }
        }
        NumTri = cont_t;  // Asigna el número total de triángulos.
        cout << "Numero de Triangulos  " << NumTri<< endl;

        //Inicializacion de las matrices
        mD.Init(NumTri, NumTri);
        mTV.Init(NumTri, NumTri);;
        mAS.Init(NumTri, NumTri);
        mPE.Init(NumTri, NumTri);

        // Inicialización de las matrices Receptor
        mDR.Init(NumReceptores, NumTri);
        mTVR.Init(NumReceptores, NumTri);
        mASR.Init(NumReceptores, NumTri);
        mPER.Init(NumReceptores, NumTri);

        int cont = 0;
        double* suma_angulos_solidos = new double[NumTri](); //arreglo para optener la sunma de las areas de los angulos solidos
        //Ciclo para el calculo de las matriz de distacia, tiempo de vuelo y angulos solidos
        for (int i = 0; i < r.NP; i++) {                      // Recorre los planos de la sala.
            for (int j = 0; j < r.p[i].NT; j++) {             // Recorre los triángulos del plano.
                int idTri1 = r.p[i].t[j].ID;                  //Obtiene el id del triangulo 1
                for (int k = 0; k < r.NP; k++) {              //Se reccorre de nuevo los planos de la sala.
                    for (int l = 0; l < r.p[k].NT; l++) {     //Recorre los triángulos del plano.
                        int idTri2 = r.p[k].t[l].ID;          //Obtiene el id del triangulo 2
                        if (i != k) {                    //verifica que ambos triángulos no pertenezcan al mismo plano.
                            mD.A[idTri1][idTri2] = r.p[i].t[j].bc.distancia(r.p[k].t[l].bc); //Calculo de la distancia de baricentro de los triangulos de la sala hacia otro triangulo de la sala NO coplanar
                            mTV.A[idTri1][idTri2] = int(1000 * mD.A[idTri1][idTri2] / V_SON); // Calculo de la distancia de baricentro de los triangulos de la sala hacia otro triangulo de la sala NO coplanar
                            mAS.A[idTri1][idTri2] = r.p[k].t[l].AnguloSolido(r.p[i].t[j].bc); // Calculo de la distancia de baricentro de los triangulos de la sala hacia otro triangulo de la sala NO coplanar
                            suma_angulos_solidos[cont] += mAS.A[idTri1][idTri2];
                        }
                    }
                }
                cont++;
            }
        }
        //Calculo de la matriz porcentaje de energia
        for (int i = 0; i < NumTri; i++) {
            for (int j = 0; j < NumTri; j++) {
                mPE.A[i][j] = mAS.A[i][j] / suma_angulos_solidos[i];
            }
        }

        cout << "Matriz distancia entre Baricentros de los triangulos" << endl;
        for (int i = 0; i < NumTri; i++) {
            for (int j = 0; j < NumTri; j++) {
                cout << fixed << setprecision(2) << mD.A[i][j] << " ";
            }
            cout << endl;
        }

        cout << "Matriz tiempo de vuelo entre Baricentros de los triangulos" << endl;
        for (int i = 0; i < NumTri; i++) {
            for (int j = 0; j < NumTri; j++) {
                cout << mTV.A[i][j] << " "; // Use setw for formatting integers
            }
            cout << endl;
        }

        cout << "Matriz porcentaje de energia de los triangulos" << endl;
        for (int i = 0; i < NumTri; i++) {
            for (int j = 0; j < NumTri; j++) {
                cout << fixed << setprecision(3) << mPE.A[i][j] << " ";
            }
            cout << endl;
        }
        //libera la memoria
        delete[] suma_angulos_solidos;


        double* suma_angulos_solidos_receptor = new double[r.NR]();

        // Ciclo para el cálculo de las matrices de distancia, tiempo de vuelo y ángulos sólidos de los receptores
        for (int i = 0; i < r.NP; i++) {                      // Recorre los planos de la sala.
            for (int j = 0; j < r.p[i].NT; j++) {             // Recorre los triángulos del plano.
                int idTri = r.p[i].t[j].ID;
                // Cálculo de las distancias, tiempo de vuelo y ángulos sólidos de los receptores
                for (int m = 0; m < r.NR; m++) {
                    mDR.A[m][idTri] = r.r[m].p.distancia(r.p[i].t[j].bc); // Cálculo de la distancia de la posición del receptor al baricentro de los triángulos de la sala
                    mTVR.A[m][idTri] = int(1000 * mDR.A[m][idTri] / V_SON); // Cálculo del tiempo de vuelo entre el receptor y el baricentro de los triángulos de la sala
                    mASR.A[m][idTri] = r.r[m].CalcularAreaDiscoProyectado(r.p[i].t[j].bc); // Cálculo del ángulo sólido entre el receptor y el baricentro de los triángulos de la sala
                    suma_angulos_solidos_receptor[m] += mASR.A[m][idTri];
                }
            }
        }

        // Cálculo de la matriz porcentaje de energía del receptor
        for (int i = 0; i < NumReceptores; i++) {
            for (int j = 0; j < NumTri; j++) {
                mPER.A[i][j] = mASR.A[i][j] / suma_angulos_solidos_receptor[i];
            }
        }

        cout << "\n\nMatriz distancia entre Baricentros de los receptores\n\n" << endl;
        for (int i = 0; i < r.NR; i++) {
            for (int j = 0; j < NumTri; j++) {
                cout << fixed << setprecision(2) << mDR.A[i][j] << " ";
            }
            cout << endl;
        }

        cout << "\n\nMatriz tiempo de vuelo entre Baricentros de los receptores\n\n" << endl;
        for (int i = 0; i < r.NR; i++) {
            for (int j = 0; j < NumTri; j++) {
                cout << mTVR.A[i][j] << " "; // Use setw for formatting integers
            }
            cout << endl;
        }

        cout << "\n\nMatriz porcentaje de energia de los receptores\n\n" << endl;
        for (int i = 0; i < r.NR; i++) {
            for (int j = 0; j < NumTri; j++) {
                cout << fixed << setprecision(3) << mPER.A[i][j] << " ";
            }
            cout << endl;
        }
        salaIniciada = true; // Indica que la habitación ha sido cargada exitosamente.
    }
}

void calcularEnergia(double x, double y, double z)
{
    if (salaIniciada)
    {
        int t_vuelo = 0; // Tiempo de vuelo del rayo
        // Coordenadas de la fuente
        o.x = x;
        o.y = y;
        o.z = z;

        s.eF = energiaFuente; //Energia de la fuente

        s.createRays(N_RAYOS);  //Numero de rayos

        double eneRayo = s.eF / s.NRAYS; // Energía Inicial de cada rayo

        cout << "\n\nEnergia inicial de cada rayo: " << eneRayo << "\n";

        reflexiones = r.RayTracing(o, s.Rays, s.NRAYS);

        mTE.Init(N_RAYOS, NUM_REBOTES); //inicializacion de la matriz de ID de triangulos en el que rebota el rayo

        mE.Init(NumTri, tM); //inicializacion de la matriz de energia

        mER.Init(r.NR, tM); //inicializacion de la matriz de energia del receptor

        // Difusión de la energía difusa en los triángulos
        for (int rayo = 0; rayo < s.NRAYS; rayo++) {
            double eneResidual = eneRayo;
            double distAcum = 0;
            for (int re = 0; re < reflexiones[rayo].N; re++) {
                int tri = reflexiones[rayo].idTriangle[re];
                distAcum += reflexiones[rayo].d[re];
                int tim = int(1000 * distAcum / V_SON);
                // Error
                //for (int j = 0; j < r.NR; j++) {
                  //  mER.A[numre][tim] = eneResidual; //Energia  del rayo al receptor
                //}
                mE.A[tri][tim] += (eneResidual * (1 - alfa) * delta); //Energia difusa en los triangulos
                eneResidual = eneResidual * (1 - alfa) * (1 - delta); //Energia incidente de los rayos
                mTE.A[rayo][re] = tri;
            }
        }

        //Transicion de energia en la matriz esapcio tiempo mE
        for (int t = 0; t < tM; t++) {
            for (int e = 0; e < NumTri; e++) {// Triángulo 1
                for (int ed = 0; ed < NumTri; ed++) { // Triángulo 2
                    //Energia de la los triangulos de la sala
                    if (e != ed) { //No estan en el mismo plano
                        //Es el instante de tiempo de la simulacion + el instate de  tiempo en que ocurre la transmision de energia en ese triangulo a hacia el siguiente triangulo
                        t_vuelo = mTV.A[e][ed] + t;
                        if (t_vuelo <= tM) {
                            mE.A[ed][t_vuelo] += (mE.A[e][t] * mPE.A[e][ed]) * (1 - alfa);
                        }
                    }
                }
                //Energía del Receptor
                for (int k = 0; k < r.NR; k++) {
                    t_vuelo = mTVR.A[k][e] + t;
                    if (t_vuelo < tM) {
                        mER.A[k][t_vuelo] += (mE.A[e][t] * mPER.A[k][e]);
                    }
                };
            }
        }

//-------------------------------Guardar en Archivo CSV----------------------------------
        // Nombre del archivo CSV
        string nombreArchivoE = "matrizTransicionEnergiaReceptores.csv";
        // Llamar a la función para guardar la matriz en el archivo CSV
        guardarMatrizCSV(mER, nombreArchivoE);
        cout << "\nMatriz guardada en " << nombreArchivoE << endl;
    }
}


void guardarMatrizCSV(MatrizResultados matriz, const string& nombreArchivo)
{
    ofstream archivo(nombreArchivo);

    if (!archivo.is_open())
    {
        cerr << "Error al abrir el archivo." << endl;
        return;
    }

    // Escribir los datos de la matriz en el archivo CSV
    for (int i = 0; i < NumReceptores; i++)
    {
        for (int j = 0; j < tM; j++)
        {
            archivo << matriz.A[i][j] << ";";
        }
        archivo << "\n";  // Nueva línea después de cada fila
    }
    archivo.close();
}
