#ifndef CHECKPOINTINPUTSTREAM_H
#define CHECKPOINTINPUTSTREAM_H

#include "StateRestore.h"
#include <fstream>
using std::istream;
using std::ios;
using namespace ProtoMol;

#ifdef HAVE_BOOST
#include <cstdlib>
#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ using ::rand; }
#endif

// the following is required to be sure the "EXPORT" works if it is used
#define CUSTOM_ARCHIVE_TYPES portable_binary_oarchive,portable_binary_iarchive
#include <boost/archive/basic_binary_iarchive.hpp>
#include <boost/archive/basic_binary_iprimitive.hpp>
#include <boost/archive/impl/basic_binary_iarchive.ipp>
#include <boost/archive/impl/basic_binary_iprimitive.ipp>
#include <boost/archive/detail/common_iarchive.hpp>
#include <boost/archive/detail/iserializer.hpp>
#include <boost/archive/binary_iarchive.hpp>
class portable_binary_iarchive :
    // don't derive from binary_iarchive !!!
    public boost::archive::binary_iarchive_impl<portable_binary_iarchive>
{
public:
  ~portable_binary_iarchive() {}
private:
  typedef portable_binary_iarchive derived_t;
  friend class boost::archive::detail::common_iarchive<derived_t>;
  friend class boost::archive::basic_binary_iarchive<derived_t>;
  friend class boost::archive::basic_binary_iprimitive<derived_t, std::istream>;
  friend class boost::archive::load_access;
  using boost::archive::binary_iarchive_impl<derived_t>::load;
 public:
  portable_binary_iarchive(std::istream & is, unsigned int flags) :
    boost::archive::binary_iarchive_impl<portable_binary_iarchive>(is, flags) {}
};
#endif


class CheckpointInputStream : 
#ifdef HAVE_BOOST
public portable_binary_iarchive
#else
public std::ifstream
#endif
{
 public:
#ifdef HAVE_BOOST
  CheckpointInputStream(std::ifstream& is, unsigned int flags)
    : portable_binary_iarchive(is, flags) {}
#else
  CheckpointInputStream(const char* filename, std::ios_base::openmode mode = std::ios_base::in)
    : std::ifstream(filename, mode) {}
#endif

  void readArchive(StateRestore& rst);
};


#endif
