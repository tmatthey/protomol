/*  -*- c++ -*-  */
#ifndef VMDSOCK_H
#define IVMDSOCK_H

#if defined(VMDSOCKINTERNAL)

#ifndef WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <sys/file.h>
#endif

namespace ProtoMol {

  namespace IMD {
    typedef struct {
      struct sockaddr_in addr; //< address of socket provided by bind() 
      int addrlen;             //< size of the addr struct 
      int sd;                  //< socket file descriptor 
    } vmdsocket;
  }
}
#endif /* VMDSOCKINTERNAL */

namespace ProtoMol {

  namespace IMD {
    int   vmdsock_init(void);
    void *vmdsock_create(void);
    int   vmdsock_bind(void *, int);
    int   vmdsock_listen(void *);
    void *vmdsock_accept(void *);  //< return new socket 
    int   vmdsock_connect(void *, const char *, int);
    int   vmdsock_write(void *, const void *, int);
    int   vmdsock_read(void *, void *, int);
    int   vmdsock_selread(void *, int);
    int   vmdsock_selwrite(void *, int);
    void  vmdsock_destroy(void *);

  }
}
#endif
