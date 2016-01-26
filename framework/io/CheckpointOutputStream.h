#ifndef CHECKPOINTOUTPUTSTREAM_H
#define CHECKPOINTOUTPUTSTREAM_H

#include "StateRestore.h"
#include "Vector3DBlock.h"
#include "ScalarStructure.h"

#include <fstream>
using std::ofstream;

#ifdef HAVE_BOOST
#include <cstdlib>
#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ using ::rand; }
#endif

// the following is required to be sure the "EXPORT" works if it is used
#define CUSTOM_ARCHIVE_TYPES portable_binary_oarchive,portable_binary_iarchive
#include <boost/archive/basic_binary_oarchive.hpp>
#include <boost/archive/basic_binary_oprimitive.hpp>
#include <boost/archive/impl/basic_binary_oarchive.ipp>
#include <boost/archive/impl/basic_binary_oprimitive.ipp>
#include <boost/archive/detail/common_oarchive.hpp>
#include <boost/archive/detail/oserializer.hpp>
#include <boost/archive/binary_oarchive.hpp>
class portable_binary_oarchive :
    // don't derive from binary_oarchive !!!
    public boost::archive::binary_oarchive_impl<portable_binary_oarchive>
{
  typedef portable_binary_oarchive derived_t;
  friend class boost::archive::detail::common_oarchive<derived_t>;
  friend class boost::archive::basic_binary_oarchive<derived_t>;
  friend class boost::archive::basic_binary_oprimitive<derived_t, std::ostream>;
  friend class boost::archive::save_access;
  using boost::archive::binary_oarchive_impl<derived_t>::save;
 public:
  portable_binary_oarchive(std::ostream & os, unsigned int flags) :
    boost::archive::binary_oarchive_impl<portable_binary_oarchive>(os, flags) {}
};
#endif

namespace ProtoMol {

class GenericTopology;
class Integrator;

class CheckpointOutputStream : 
#ifdef HAVE_BOOST
public portable_binary_oarchive
#else
public std::ofstream
#endif
{
 public:
#ifdef HAVE_BOOST
  CheckpointOutputStream(std::ostream& os, unsigned int flags)
    : portable_binary_oarchive(os, flags) {}
#else
  CheckpointOutputStream(const char* filename, std::ios_base::openmode mode = std::ios_base::out)
    : std::ofstream(filename, mode) {}
#endif


  void writeArchive(GenericTopology* topo,
		    Vector3DBlock positions,
		    Vector3DBlock velocities,
		    Integrator* integrator,
		    ScalarStructure scalar);
};
}

#endif
