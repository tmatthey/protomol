/*  -*- c++ -*-  */
#ifndef WRITER_H
#define WRITER_H

#include "File.h"

namespace ProtoMol {

  //_________________________________________________________________ Writer
  /**
   * Base clas of writers
   */
  class Writer : public File {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    Writer();
    explicit Writer(const std::string& filename);
    /// To open with special file flags, std::ios::out is set
    explicit Writer(std::ios::openmode mode);
    /// To open with special file flags, std::ios::out is set
    Writer(std::ios::openmode mode, const std::string& filename);
  public:
    virtual ~Writer();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Writer
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /**
     * A 1-line comment, which might added to the output by the actual 
     * writer implementation; no new line!
     */
    void setComment(const std::string& comment="");

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:

  };

  //______________________________________________________________________ INLINES
}

#endif /* WRITER_H */
