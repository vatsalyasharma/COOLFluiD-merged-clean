#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRhsJacobCellParalution_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRhsJacobCellParalution_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeCUDA/FVMCC_ComputeSourceRhsJacobCell.hh"
#include "FiniteVolume/KernelData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class CellConn;
  }
  
  namespace Numerics {

    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS using
 * standard cell center FVM schemes with CUDA bindings
 *
 * @author Andrea Lani
 *
 */

//This template include one more parameter, the source term
template <typename SCHEME, typename PHYSICS, typename SOURCE, typename POLYREC, typename LIMITER, CFuint NB_BLOCK_THREADS>
class FVMCC_ComputeSourceRhsJacobCellParalution: public FVMCC_ComputeSourceRhsJacobCell<SCHEME, PHYSICS, SOURCE, POLYREC, LIMITER, NB_BLOCK_THREADS> {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeSourceRhsJacobCellParalution(const std::string& name) :
    FVMCC_ComputeSourceRhsJacobCell<SCHEME, PHYSICS, SOURCE, POLYREC, LIMITER, NB_BLOCK_THREADS>(name) {}
  
  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeSourceRhsJacobCellParalution() {} 
 
  /**
   * Execute Processing actions
   */
  virtual void execute();
  
}; // class FVMCC_ComputeSourceRhsJacobCellParalution

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeSourceRhsJacobCellParalution_hh
