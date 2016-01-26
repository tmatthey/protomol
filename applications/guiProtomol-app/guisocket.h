#ifndef __GUISOCKET_H__
#define __GUISOCKET_H__

#define COMM_MAGIC   0xF01DBAAD
#define COMM_NNCOORD 050363

#define COMM_VERSION 1
#define COMM_PORT    0xCE11

#ifdef _MSC_VER

typedef __int8            int8_t;
typedef __int16           int16_t;
typedef __int32           int32_t;
typedef __int64           int64_t;
typedef unsigned __int8   uint8_t;
typedef unsigned __int16  uint16_t;
typedef unsigned __int32  uint32_t;
typedef unsigned __int64  uint64_t;

#else
#include <stdint.h>
#endif

typedef struct {uint64_t x[2];} uint128_t;

typedef struct {
  uint8_t type[4];
  float charge;
  float radius;
} FAH_ATOM;

typedef struct {
  uint32_t a; /* rule: a < b */
  uint32_t b;
} FAH_BOND;

typedef struct {
  float x;
  float y;
  float z;
} FAH_XYZ;

typedef struct {
  uint32_t magic;
  uint32_t version;
  uint8_t  name[64];
  int64_t  timestamp;
  uint64_t iterations;
  uint32_t frames;
  uint32_t atom_count;
  uint32_t bond_count;
} FAH_INFO;

typedef struct {
  uint32_t magic;
  uint32_t version;
  int64_t  timestamp;
  uint64_t iterations_done;
  uint32_t frames_done;
  float    energy;
  float    temperature;
} FAH_CURRENT;

typedef enum {
  GS_NO_REQUEST,
  GS_META_REQUEST,
  GS_COORD_REQUEST,
} gs_request_t;

typedef struct {
  FAH_INFO info;
  FAH_CURRENT current;
  FAH_ATOM *atoms;
  FAH_BOND *bonds;
  FAH_XYZ *xyz;
} gs_data_t;

#if defined WIN32
#define NOMINMAX
#include <winsock.h>

#else
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#include <netdb.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#endif

#include <time.h>

static gs_request_t gs_request = GS_NO_REQUEST;
static gs_data_t gs_data;
static int gs_shutdown = 0;
static unsigned short com_port;
static unsigned int com_port_range;
static uint32_t no_new_coord = COMM_NNCOORD;

// Redefine some types and constants based on OS
#ifdef WIN32
typedef int socklen_t;           // Unix socket length
#else
typedef int SOCKET;
#define INVALID_SOCKET -1        // WinSock invalid socket
#define SOCKET_ERROR   -1        // Basic WinSock error
#define closesocket(s) close(s); // Unix uses file descriptors, WinSock doesn't
#endif

/* This code sets up a server. When the server gets a connection, 
 * it reads a 4 byte request code. We then set a flag with the value of 
 * the requested data. This flag is examined periodically by the 
 * core and the core may choose to fill in values when it is good
 * and ready. When it does so, the flag gets reset indicating
 * that the data has been filled in. If the client can not wait,
 * then the core wrapper has no choice but to return old values
 * since the science code may be busy. For now, I'm assuming the 
 * client code can smoothly handle delays from the core. So we 
 * just use one mutex that indicates a request has been made. 
 * While this mutex is locked, the core wrapper code will not 
 * touch the coordinates etc.
 */


#ifdef WIN32
static HANDLE gs_thread;
static HANDLE gs_mutex;
#else
static pthread_t gs_thread;
static pthread_mutex_t gs_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif


static void gs_lock(void) {
#ifdef WIN32
  WaitForSingleObject(gs_mutex, INFINITE);
#else
  pthread_mutex_lock(&gs_mutex);
#endif
}

static void gs_unlock(void) {
#ifdef WIN32
  ReleaseSemaphore(gs_mutex, 1, 0);
#else
  pthread_mutex_unlock(&gs_mutex);
#endif
}


/* This will be called once at startup. The science code tells
 * us how big the system is and gives us a callback function 
 * that we can use for filling in the data structures
 * 
 * The science code will give us the number of atoms and bonds
 * Note that this gets called _after_ the science thread has been 
 * launched
 */
static void gs_init_data(int natoms, int nbonds,  char *name, int frames) {
  /* Fill in stuff that the science code need not bother with
   * Also, init to non-junk values, since the viewer may update
   * before the science code has had a chance to fill in anything
   */
  gs_data.info.version = COMM_VERSION;
  gs_data.info.magic   = COMM_MAGIC;

  strncpy((char *)gs_data.info.name, name, 64);
  gs_data.info.name[63] = 0;

  gs_data.info.timestamp = time(0);

  gs_data.info.frames     = frames;
  gs_data.info.iterations = 0;

  gs_data.info.atom_count = natoms;
  gs_data.info.bond_count = nbonds;
    
  gs_data.current.energy      = 0.0f;
  gs_data.current.temperature = 0.0f;
  gs_data.current.magic   = COMM_MAGIC;
  gs_data.current.version = COMM_VERSION;
  gs_data.current.iterations_done = 0;
  gs_data.current.frames_done = 0;
  gs_data.current.timestamp = gs_data.info.timestamp;

  gs_data.atoms = new FAH_ATOM[natoms];
  gs_data.bonds = new FAH_BOND[nbonds];
  gs_data.xyz   = new FAH_XYZ[natoms];
}

static void gs_free_data(void) {
  delete [] gs_data.atoms;
  delete [] gs_data.bonds;
  delete [] gs_data.xyz;

  gs_data.atoms = 0;
  gs_data.bonds = 0;
  gs_data.xyz = 0;
}

static SOCKET gs_send(SOCKET socket, char *data, unsigned length) {
  if (send(socket, data, length, 0) != (int)length) {
    //printf("Network error sending data\n");
    return -1;
  }

  return 0;
}

static int gs_send_metadata(SOCKET socket) {
  return
    gs_send(socket, (char *)&gs_data.info, sizeof(FAH_INFO)) || //100
    gs_send(socket, (char *)gs_data.atoms,
            sizeof(FAH_ATOM) * gs_data.info.atom_count) ||
    gs_send(socket, (char *)gs_data.bonds, sizeof(FAH_BOND) * gs_data.info.bond_count);
}

static int gs_send_coordinates(SOCKET socket) {
  return
    gs_send(socket, (char *)&gs_data.current, sizeof(FAH_CURRENT)) || //36
    gs_send(socket, (char *)gs_data.xyz, sizeof(FAH_XYZ) * gs_data.info.atom_count);
}

static int gs_send_no_new_coord(SOCKET socket) {
  return
    gs_send(socket, (char *)&no_new_coord, sizeof(uint32_t));
}

#ifdef WIN32
static DWORD WINAPI gs_server(LPVOID param)
#else
static void *gs_server(void *param)
#endif
{
    struct sockaddr_in addr, r_addr;
    SOCKET lsocket = INVALID_SOCKET;
    SOCKET asocket = INVALID_SOCKET;
    socklen_t len = sizeof(r_addr);   // The length of our remote address
      
    if ((lsocket = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) == INVALID_SOCKET) {
        //printf("Invalid socket!\n");
        return 0;
    }
    
  // Setup the local address variable
  memset((void *)&addr, 0, sizeof(addr));
  addr.sin_family      = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_ANY);
  addr.sin_port        = htons(com_port);//COMM_PORT);
  
  // Name the local socket
  if (bind(lsocket, (struct sockaddr *)&addr, sizeof(addr)) == SOCKET_ERROR) {
    //try one port up if failed (use for dual/quad core machines)
    if(com_port_range > 1){
        int bret = SOCKET_ERROR;
        for(unsigned int i=1;i<com_port_range && bret == SOCKET_ERROR;i++){
            addr.sin_port = htons(com_port+i);
            bret = bind(lsocket, (struct sockaddr *)&addr, sizeof(addr));
        }
        if(bret == SOCKET_ERROR){
            perror("Could not bind to local socket");
            closesocket(lsocket);
            return 0;
        }
    }else{
        perror("Could not bind to local socket");
        closesocket(lsocket);
        return 0;
    }
    //addr.sin_port = htons(com_port+1);
    //if (bind(lsocket, (struct sockaddr *)&addr, sizeof(addr)) == SOCKET_ERROR) {
    //    perror("Could not bind to local socket");
    //    closesocket(lsocket);
    //    return 0;
    //}
  }
  
  // Set the socket to listen for a connection
  if (listen(lsocket, SOMAXCONN) == SOCKET_ERROR) {
    perror("Could not listen to local socket");
    closesocket(lsocket);
    return 0;
  }

  // Wait for a connection
// While not close connection
  while (!gs_shutdown)
    //printf("GUI Server main loop\n");
    if ((asocket = accept(lsocket, (struct sockaddr *)&r_addr, &len)) !=
        INVALID_SOCKET) {

      //printf("GUI Server connected\n");

      int close_connection = 0;
// While not close connection
      while (!close_connection && !gs_shutdown) {
        unsigned request;
        //printf("GUI Server connection open\n");
// Recieve data
        if (recv(asocket, (char *)&request, 4, 0) == 4) {

// Switch on request
         switch (request) {
          case 0:
          case 1:

            //gs_lock();
            //if (request == 0) gs_request = GS_META_REQUEST;
            //else gs_request = GS_COORD_REQUEST;
            //gs_unlock();

            //printf("GUI Server data requested %i \n",request);
            if (gs_shutdown) break;

            //gs_lock();
            //printf("GUI Server sending data\n");
// Send data
            if (request == 0) {
                if (gs_send_metadata(asocket)){
                    close_connection = 1;
                    //printf("GUI Server meta close=1\n");
                }
            } else{
                gs_lock();
                if(gs_request == GS_COORD_REQUEST){ //no new data?
                    if(gs_send_no_new_coord(asocket)){
                        close_connection = 1;
                    }
                }else{
                    gs_request = GS_COORD_REQUEST;
                    if (gs_send_coordinates(asocket)){
                        close_connection = 1;
                        //printf("GUI Server data close=1\n");
                    }
                }
                gs_unlock();
            }
// end data sent
            //gs_unlock();
            break;

          case 99: close_connection = 1; break;

          default:
            //printf("GUI Server: bad request %u\n", request);
            close_connection = 1;
         }
// end switch on request
        }else{
            close_connection = 1;
            //printf("GUI Server end session\n");
        }
// end recieve data
      }
// While not close connection
      closesocket(asocket);
    }
// end while not close connection
  
  closesocket(lsocket);

  return 0;
}

static int gs_start_server(unsigned short c_port, unsigned int c_p_range) {
    com_port = c_port;	//save port
    com_port_range = c_p_range;
#ifdef WIN32
  WSADATA wsa;
  int err;

  if ((err = WSAStartup(MAKEWORD(2, 2), &wsa)) != NO_ERROR) {
    //printf("Error starting winsock\n");
    WSACleanup();
    return -1;
  }

  if (LOBYTE(wsa.wVersion) != 2 || HIBYTE(wsa.wVersion) != 2) {
    //printf("Error need winsock version 2.2\n");
    WSACleanup();
    return -1;
  }
#endif

  gs_shutdown = 0;

#ifdef WIN32
  DWORD tid;

  gs_mutex = CreateSemaphore(0, 1, 1, 0);

  if (!(gs_thread = CreateThread(0, 0, gs_server, 0, 0, &tid))) {
    //printf("Error starting GUI server thread.\n");
    return -1;
  }

#else

  if (pthread_create(&gs_thread, 0, gs_server, 0)) {
    printf("Error starting GUI server thread.\n");
    return -1;
  }
#endif

  return 0;
}

gs_data_t *gs_init(int natoms, int nbonds, char *name, int frames) {
  gs_init_data(natoms, nbonds, name, frames);

  if (gs_start_server(0xCE11,1) == -1) return 0;

  return &gs_data;
}

gs_request_t gs_get_request(void) {
  return gs_request;
}

void gs_start_update(void) {
  gs_lock();
}

void gs_end_update(void) {
  //printf("GUI data served\n");
  gs_request = GS_NO_REQUEST;
  gs_unlock();
}

void gs_cleanup(void) {
  gs_shutdown = 1;

  // Wait for server thread
#ifdef WIN32
  if (gs_thread)
    WaitForSingleObject(gs_thread, 5000); // Upto 5sec
#endif

  gs_free_data();

#ifdef WIN32
  if (gs_thread) CloseHandle(gs_thread);
  gs_thread = 0;
  if (gs_mutex) CloseHandle(gs_mutex);
  gs_mutex = 0;

  WSACleanup();
#endif
}

static int gui_meta( GenericTopology *top ) {
    unsigned int i;	
    int element;
    float mass;
    float myRadius;

   //report << plain << " gui_bonds" <<endl;
    
   //bonds
    for ( i = 0; i < top->bonds.size(); i++ ) {
        //FAH requires a < b
        if (top->bonds[i].atom1 < top->bonds[i].atom2)
        {
            gs_data.bonds[i].a = top->bonds[i].atom1;
            gs_data.bonds[i].b = top->bonds[i].atom2;
        }else
        {
            gs_data.bonds[i].a = top->bonds[i].atom2;
            gs_data.bonds[i].b = top->bonds[i].atom1;
        }
    }
    //atoms
    //report << plain << " gui_atoms" <<endl;
    for ( i = 0; i < (top->atoms).size(); i++ ) {
        //Determine element by mass
        //This will not work with united atom models or with heavy hydrogens etc
        element = 0; //unknown
        myRadius = 2.5f;
        gs_data.atoms[i].type[0] = 'Q';
        mass = top->atoms[i].scaledMass;
        if ( mass < 1.2 && mass >= 1.0 ) //hydrogen
        {
            element = 1;
            gs_data.atoms[i].type[0] = 'H';
            myRadius = 1.2f;
        }
        else if ( mass > 11.8 && mass < 12.2 ) //carbon
        {
            element = 6;
            gs_data.atoms[i].type[0] = 'C';
            myRadius = 1.7f;
        }
        else if ( mass > 14.0 && mass < 15 )   //nitrogen
        {
            element = 7;
            gs_data.atoms[i].type[0] = 'N';
            myRadius = 1.55f;
        }
        else if ( mass > 15.5 && mass < 16.5 ) //oxygen
        {
            element = 8;
            gs_data.atoms[i].type[0] = 'O';
            myRadius = 1.52f;
        }
        else if ( mass > 31.5 && mass < 32.5 ) //sulphur
        {
            element = 16;
            gs_data.atoms[i].type[0] = 'S';
            myRadius = 1.85f;
        }
        else if ( mass > 29.5 && mass < 30.5 ) //phosphorus
        {
            element = 15;
            gs_data.atoms[i].type[0] = 'P';
            myRadius = 1.9f;
        }

        gs_data.atoms[i].charge = top->atoms[i].scaledCharge;
        gs_data.atoms[i].radius = myRadius * 0.5;

    }

    return 1;
}

static int gui_coords( Vector3DBlock *pos, int step, int mode ) {

    //report << plain << " gui_coords" <<endl;
    Real x, y, z, sz; 

    x = y = z = 0.0;
    sz = pos->size();
    for (unsigned int i = 0; i < pos->size(); i++ ) {
        x += (*pos)[i].x;
        y += (*pos)[i].y;
        z += (*pos)[i].z;
    }
    x /= sz; y /= sz; z /= sz;
    for (unsigned int i = 0; i < pos->size(); i++ ) {
        gs_data.xyz[i].x = (float)((*pos)[i].x - x);
        gs_data.xyz[i].y = (float)((*pos)[i].y - y);
        gs_data.xyz[i].z = (float)((*pos)[i].z - z);
    }
    //progress
    gs_data.current.frames_done = step;
    gs_data.current.timestamp = time(0);
    gs_data.current.iterations_done = mode;
    //
    return 1;
}


#endif
