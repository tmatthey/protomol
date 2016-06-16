#include "IMDElf.h"
#include "Timer.h"
#include "Report.h"
#include "vmdsock.h"

using namespace ProtoMol::Report;
using namespace ProtoMol::IMD;

namespace ProtoMol
{
	//_________________________________________________________________ IMDElf
	IMDElf::IMDElf(): trate(1), timeout(1000), wait_imd(0), sock(NULL), clientsock(NULL)
	{
	}

	IMDElf::~IMDElf()
	{
		if (sock)
			vmdsock_destroy(sock);
		if (clientsock)
			vmdsock_destroy(clientsock);
	}

	int IMDElf::listen(int defport)
	{
		sock = vmdsock_create();
		port = find_free_port(sock, defport);

		if (port < 0)
		{
			vmdsock_destroy(sock);
			return 1;
		}
		vmdsock_listen(sock);
		return 0;
	}

	int IMDElf::hookup()
	{
		clientsock = 0;

		if (!clientsock)
		{
			if (wait_imd)
			{
				report << plain << "Waiting for IMD connection..." << endr;
				do
				{
					rc = vmdsock_selread(sock, 3600);
				}
				while (rc <= 0);
			}
			else
			{
				rc = vmdsock_selread(sock, 0);
			}
			//while ((rc = vmdsock_selread(sock, 3600)) <= 0);

			if (rc > 0)
			{
				clientsock = vmdsock_accept(sock);
				if (!clientsock)
				{
					report << recoverable << "IMDElf: Accept failed" << endr;
					return 1;
				}
				else
				{
					if (!connect(clientsock))
					{
						report << recoverable << "IMDElf: IMD handshake failed" << endr;
						vmdsock_destroy(clientsock);
						clientsock = NULL;
						return 1;
					}
				}
			}
		}
		return 0;
	}

	int IMDElf::send_coords(const Vector3DBlock& c)
	{
		//First convert the coordiantes to a giant array of floats so that VMD
		//can take it.
		float* coords = new float[c.size() * 3];
		//float coords[(const int) c.size() * 3];

		for (unsigned int i = 0; i < c.size(); ++i)
		{
			coords[i * 3] = c[i].x;
			coords[i * 3 + 1] = c[i].y;
			coords[i * 3 + 2] = c[i].z;
		}

		int rc = imd_send_fcoords(clientsock, c.size() * 3, (const float*) coords);
		delete [] coords;

		return rc;
	}

	int IMDElf::get_haptics(Vector3DBlock& h)
	{
		int length;
		bool kill_switch = false;

		while (!kill_switch && (vmdsock_selread(clientsock, 0) > 0))
		{ // Drain the socket
			int type = imd_recv_header(clientsock, &length);

			switch (type)
			{
			case IMD_MDCOMM:
				vmd_atoms = new int32[length];
				vmd_forces = new float[3 * length];
				h.resize(length); //Resize the haptic forces we will return
				if (imd_recv_mdcomm(clientsock, length, vmd_atoms, vmd_forces))
				{
					report << recoverable << "Error reading MDComm forces" << endr;
					return -1;
				}

				for (int i = 0; i < length; ++i)
				{
					h[vmd_atoms[i]].x = vmd_forces[i * 3];
					h[vmd_atoms[i]].y = vmd_forces[i * 3 + 1];
					h[vmd_atoms[i]].z = vmd_forces[i * 3 + 2];
				}

				delete [] vmd_atoms;
				delete [] vmd_forces;
				break;
			case IMD_TRATE:
				if (length > 0)
				{
					report << plain << "Setting transfer rate to client at " << length << endr;
					trate = length;
				}
				else //Invalid transfer rate
					report << recoverable << "Invalid transfer rate.  Rate unchanged." << endr;
				break;
			case IMD_DISCONNECT:
				report << plain << "Remote client has detached from simulation." << endr;
				vmdsock_destroy(clientsock);
				clientsock = 0;
				kill_switch = true;
				break;
			case IMD_ENERGIES:
				IMDEnergies junk;
				imd_recv_energies(clientsock, &junk);
				break;
			case IMD_FCOORDS:
				vmd_forces = new float[3 * length];
				imd_recv_fcoords(clientsock, length, vmd_forces);
				delete [] vmd_forces;
				break;
			case IMD_IOERROR:
				vmdsock_destroy(clientsock);
				clientsock = NULL;
				report << recoverable << "IMD connection lost" << endr;
				break;
			default:
				report << recoverable << "Remote client sent invalid header" << endr;
			}
		}

		return 0;
	}

	bool IMDElf::client_alive()
	{
		if (clientsock == 0)
			return false;
		return true;
	}

	int IMDElf::find_free_port(void* sock, int defport)
	{
		if (vmdsock_bind(sock, defport) == 0) return defport; // success
		for (int port = 1025; port < 4096; port++)
			if (vmdsock_bind(sock, port) == 0) return port;
		return -1;
	}

	int IMDElf::connect(void* s)
	{
		if (imd_handshake(s))
		{
			//iout << iWARN << "IMD handshake failed\n" << endi;
			report << recoverable << "IMD handshake failed" << endr;
			return 0;
		}

		// Wait a second, then see if VMD has responded.
		//double t = CmiWallTimer();
		//int i = 10000000;
		//while (i > 0)
		//  --i;

		Timer t;
		t.start();
		while (1000 * t.getTime().getUserTime() < timeout);

		//While (CmiWallTimer()-t < 1.0);
		int32 length;
		if (vmdsock_selread(s, 0) != 1 || imd_recv_header(s, &length) != IMD_GO)
		{
			report << recoverable << "Incompatible Interactive MD, use VMD v1.4b2 or higher" << endr;
			return 0;
		}
		return 1;
	}
}
