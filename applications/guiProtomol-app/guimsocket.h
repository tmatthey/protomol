#ifndef __GUISOCKET_H__
#define __GUISOCKET_H__

#define COMM_MAGIC   0xF01DBAAD
#define COMM_NNCOORD 050363

#define COMM_VERSION 2
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
#include <errno.h>
#include <arpa/inet.h>
#include <fcntl.h>
#endif

#include <time.h>

static gs_request_t gs_request = GS_COORD_REQUEST;//Start with request so no wait. //GS_NO_REQUEST;
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

#define NUMCONN 5
static SOCKET connectlist[NUMCONN];
bool connectlist_sent_curr[NUMCONN];
static fd_set socks; 

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
  int sendRet;
  unsigned numSent = 0;
  //non-blocking so wait for buffer
  while(numSent < length){ 
    sendRet = send(socket, &data[numSent], length-numSent, 0);
    if(sendRet < 0){
        //printf("Network error sending data\n");
#ifndef WIN32
        if(errno != EWOULDBLOCK) return -1;
#else
        if(WSAGetLastError() != WSAEWOULDBLOCK) return -1;
#endif
    }else{
        numSent += sendRet;
    }
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
  }
  
  // Set the socket to listen for a connection
  if (listen(lsocket, SOMAXCONN) == SOCKET_ERROR) {//
    perror("Could not listen to local socket");
    closesocket(lsocket);
    return 0;
  }

  // Wait for a connection
  //initialize connection list
  for(int listnum = 0; listnum < NUMCONN; listnum++){
      connectlist[listnum] = INVALID_SOCKET;
      connectlist_sent_curr[listnum] = true;	//don't send until real data
  }

  //####Main server loop###################################################################
  // While not close connection
  while (!gs_shutdown){

    int highsock=lsocket;
    int close_connection = 0;
    // Wait for a connection
    //initialize file descriptors
    FD_ZERO(&socks);
    FD_SET(lsocket,&socks);
    for(int listnum = 0; listnum < NUMCONN; listnum++){
        if(connectlist[listnum] != INVALID_SOCKET){
            FD_SET(connectlist[listnum],&socks);
            if(connectlist[listnum] > highsock)
                highsock = connectlist[listnum];
        }
    }
    //wait for socket activity
    int readsocks = select(highsock+1, &socks, (fd_set *) 0, (fd_set *) 0, NULL);//&timeout, 
    if(readsocks < 0) {
        perror("select");
        return 0;
    }
    if(readsocks){
        //activity on new connection?
        if(FD_ISSET(lsocket,&socks)){	
            int socket_error = 0;
            if((asocket = accept(lsocket, (struct sockaddr *)&r_addr, &len)) != INVALID_SOCKET){
                //set non blocking
#ifndef WIN32
                int opts;
                if((opts = fcntl(asocket,F_GETFL)) >= 0){
                    opts = (opts | O_NONBLOCK);
                    if(fcntl(asocket,F_SETFL,opts) >= 0){
#else             
                u_long on = 1;
                if(ioctlsocket(asocket,FIONBIO, &on) != SOCKET_ERROR){
#endif
                        //put in list
                        for(int listnum = 0; (listnum < NUMCONN) && (asocket != -1); listnum ++){
                            if (connectlist[listnum] == INVALID_SOCKET){
                                printf("Connection accepted:   FD=%d; Slot=%d\n", asocket, listnum);
                                connectlist[listnum] = asocket;
                                asocket = -1;
                            }
                        }
                        if (asocket != -1) {
                            // No room left in the queue
                            printf("No room left for new client.\n");
                            socket_error = 1;
                        }
#ifndef WIN32
                    }else{
                        printf("Error on fcntl(F_SETFL)\n");
                        socket_error = 1;
                    }
                }else{
                    printf("Error on fcntl(F_GETFL)\n");
                    socket_error = 1;
                }
#else             
                }else{
                    printf("Error on ioctl\n");
                    socket_error = 1;
                }
#endif
            }else{
                printf("Error on accept\n");
                socket_error = 1;
            }
            if(socket_error) closesocket(asocket);
        }
        //handle data for all entries in queue
        unsigned request;
        for(int listnum = 0; listnum < NUMCONN; listnum++){
            if(connectlist[listnum] != INVALID_SOCKET && FD_ISSET(connectlist[listnum],&socks)){
                close_connection = 0;
                if(recv(connectlist[listnum], (char *)&request, 4, 0) == 4) {
                     // Switch on request
                     switch (request) {
                          case 0:	if (gs_shutdown) break;	//send metadata
                                    if (gs_send_metadata(connectlist[listnum])){
                                        close_connection = 1;
                                        //printf("GUI Server meta close=1\n");
                                    }
                                    break;
                          case 1:   if (gs_shutdown) break;	//send data
                                    gs_lock();
                                    //Globally new daya? Required in case version 1 client in set
                                    if(gs_request != GS_COORD_REQUEST)
                                        for(int clist=0;clist < NUMCONN;clist++) connectlist_sent_curr[clist] = false;	//set data available to all clients
                                    gs_request = GS_COORD_REQUEST;
                                    if (gs_send_coordinates(connectlist[listnum])){
                                        close_connection = 3;
                                        //printf("GUI Server data close=1\n");
                                    }
                                    gs_unlock();
                                    break;
                          case 2:   if (gs_shutdown) break;	//send data if new available
                                    gs_lock();
                                    if(gs_request == GS_COORD_REQUEST && connectlist_sent_curr[listnum]){ //no new data? AND this connection has sent it
                                    //if(gs_request == GS_COORD_REQUEST){ //no new data?
                                        gs_unlock();
                                        if(gs_send_no_new_coord(connectlist[listnum])){
                                            close_connection = 2;
                                        }
                                    }else{
                                        //Globally new daya, or just this client?
                                        if(gs_request != GS_COORD_REQUEST)	//global
                                            for(int clist=0;clist < NUMCONN;clist++) connectlist_sent_curr[clist] = false;	//set data available to all clients
                                        //service request
                                        gs_request = GS_COORD_REQUEST;
                                        connectlist_sent_curr[listnum] = true;	//flag this client has had new data
                                        if (gs_send_coordinates(connectlist[listnum])){
                                            close_connection = 3;
                                            //printf("GUI Server data close=1\n");
                                        }
                                        gs_unlock();
                                    }
                                    break;

                          case 99: close_connection = 4; break;

                          default:
                            //printf("GUI Server: bad request %u\n", request);
                            close_connection = 5;
                     }
                }else{
                    close_connection = 6;
                }
                if(close_connection){
                    printf("Connection closed:   FD=%d; Slot=%d; Err=%d;\n", connectlist[listnum],listnum,close_connection);
                    closesocket(connectlist[listnum]);
                    connectlist[listnum] = INVALID_SOCKET;
                }
            }
        }
        //end handle data
    }

  }
  //####End Main server loop###############################################################
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
        // Determine diamiter/ set name
        int atomNameLen = top->atoms[i].name.length();
        //for(int j=0;j<4 && j<atomNameLen;j++) gs_data.atoms[i].type[j] = (top->atoms[i].name.c_str())[j];
        for(int j=0;j<min(atomNameLen,4);j++) gs_data.atoms[i].type[j] = (top->atoms[i].name.c_str())[j];
        if(atomNameLen < 4)
            for(int j=atomNameLen;j<4;j++) gs_data.atoms[i].type[j] = 0;
        switch(gs_data.atoms[i].type[0]){
            case 'H':	myRadius = 1.2f;
                        break;
            case 'C':	myRadius = 1.7f;
                        break;
            case 'N':	myRadius = 1.55f;
                        break;
            case 'O':	myRadius = 1.52f;
                        break;
            case 'S':	myRadius = 1.85f;
                        break;
            case 'P':	myRadius = 1.9f;
                        break;
            default:	myRadius = 1.9f;
                        break;
        }
        //add charge
        gs_data.atoms[i].charge = top->atoms[i].scaledCharge;
        gs_data.atoms[i].radius = myRadius * 0.5;

    }

    return 1;
}

static int gui_coords( Vector3DBlock *pos, int step, int mode, GenericTopology *topo ) {

    //report << plain << " gui_coords" <<endl;
    Real x, y, z, sz; 
    Vector3DBlock posMi;

    posMi.resize(pos->size());
    //
    x = y = z = 0.0;
    sz = pos->size();
    for (unsigned int i = 0; i < pos->size(); i++ ) {
        posMi[i] = topo->minimalPosition((*pos)[i]);
        x += posMi[i].x;
        y += posMi[i].y;
        z += posMi[i].z;
    }
    x /= sz; y /= sz; z /= sz;
    for (unsigned int i = 0; i < pos->size(); i++ ) {
        gs_data.xyz[i].x = (float)(posMi[i].x - x);
        gs_data.xyz[i].y = (float)(posMi[i].y - y);
        gs_data.xyz[i].z = (float)(posMi[i].z - z);
    }
    //progress
    gs_data.current.frames_done = step;
    gs_data.current.timestamp = time(0);
    gs_data.current.iterations_done = mode;
    //
    return 1;
}


#endif
