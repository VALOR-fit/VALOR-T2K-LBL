//____________________________________________________________________________
/*!

\enum      t2kana::conv::ConfLevel_t

\brief     Confidence level

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group

*/
//____________________________________________________________________________

#ifndef __valor_conf_level_h__
#define __valor_conf_level_h__

namespace t2kana {

namespace conv {

typedef enum EConfLevel {
   kCLUndefined = -1, 
   kCL68, 
   kCL90, 
   kCL95, 
   kCL99
} 
ConfLevel_t; 

} // conv

namespace util {


} // util

} // t2kana
  
#endif 

