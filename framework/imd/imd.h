/* -*- c++ -*- */
#ifndef IMD_H
#define IMD_H

#include "protomol.h"
#include "typeSelection.h"
namespace ProtoMol {

  namespace IMD {

    typedef TypeSelection::Int<4>::type   int32;

    enum IMDType {
      IMD_DISCONNECT,
      IMD_ENERGIES, 
      IMD_FCOORDS,   
      IMD_GO,
      IMD_HANDSHAKE, 
      IMD_KILL,      
      IMD_MDCOMM,    
      IMD_PAUSE,
      IMD_TRATE,
      IMD_IOERROR
    };

    typedef struct {
      int32 tstep;
      float T;
      float Etot;
      float Epot;
      float Evdw;
      float Eelec;
      float Ebond;
      float Eangle;
      float Edihe;
      float Eimpr;
    } IMDEnergies;

    /// Send simple messages - these consist of a header with no subsequent data
    int   imd_disconnect(void *);
    int   imd_pause(void *);
    int   imd_kill(void *);
    int   imd_handshake(void *);
    int   imd_trate(void *, int32);

    // Send data
    int   imd_send_mdcomm(void *, int32, const int32 *, const float *);
    int   imd_send_energies(void *, const IMDEnergies *);
    int   imd_send_fcoords(void *, int32, const float *);


    /** 
     *recv_handshake returns 0 if server and client have the same relative 
     * endianism; returns 1 if they have opposite endianism, and -1 if there
     * was an error in the handshake process.
     */
    int imd_recv_handshake(void *);

    IMDType imd_recv_header(void *, int32 *);
    int imd_recv_mdcomm(void *, int32, int32 *, float *);
    int imd_recv_energies(void *, IMDEnergies *);
    int imd_recv_fcoords(void *, int32, float *);

  }
}
#endif

