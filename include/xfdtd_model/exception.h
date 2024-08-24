#ifndef __XFDTD_MODEL_EXCEPTION_H__
#define __XFDTD_MODEL_EXCEPTION_H__

#include <xfdtd/exception/exception.h>

namespace xfdtd::model {

class XFDTDModelException : public XFDTDException {
 public:
  using XFDTDException::XFDTDException;
};

}  // namespace xfdtd::model

#endif  // __XFDTD_MODEL_EXCEPTION_H__
