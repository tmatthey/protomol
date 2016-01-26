/*  -*- c++ -*-  */
#ifndef IMDELF
#define IMDELF

#include "Vector3DBlock.h"
#include "imd.h"

namespace ProtoMol {
  /**
   * This is the IMD Elf, who is charged with maintaining all aspects of an
   * IMD connection with Protomol.  Is the Elf not nifty?  Pay homage to the
   * Elf.
   */

  //_________________________________________________________________ IMDElf

  class IMDElf {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    IMDElf();
    ~IMDElf();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Public IMDELF methods
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    int listen(int);
    int hookup();
    int send_coords(const Vector3DBlock&);
    int get_haptics(Vector3DBlock&);
    bool client_alive();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Public IMDElf data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int trate;
    int timeout;
    int wait_imd;

  private:
    void *sock;
    void *clientsock;
    int rc;
    int port;
    IMD::int32* vmd_atoms;
    float* vmd_forces;

    int find_free_port(void *, int);
    int connect(void *);  /// pass socket handle
  };
}
#endif
