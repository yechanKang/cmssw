#ifndef LINEARCALIBRATOR_HH_
#define LINEARCALIBRATOR_HH_

#include "RecoParticleFlow/PFClusterTools/interface/Calibrator.hh"

#include "TMatrixD.h"
#include "TVectorD.h"

/*
 * LinearCalibrator Class
 * 
 * This class implements the simple "linear" calibration for the "a,b,c" coefficients.
 * It extends Calibrator.
 * 
 * Calibrations are given i.t.o,
 * 		E_calib = a + b * det_1 + c * det_2 + ...
 * 
 * doOffset is set to true by default: this 'a' term accounts for threshold effects.
 */
namespace pftools {
class LinearCalibrator : public Calibrator {
public:
	LinearCalibrator();
	virtual ~LinearCalibrator();

	virtual std::map<DetectorElement*, double>
			getCalibrationCoefficients() throw(MinimiserException&);

	/*
	 * Note: covariant return type w.r.t. Calibrator class: the overloading has changed 
	 * the return type but this IS allowed with modern compilers.
	 * See documentation in Calibrator.hh
	 */
	LinearCalibrator* clone() const;
	LinearCalibrator* create() const;

protected:
	
	LinearCalibrator(const LinearCalibrator& lc);
	/*
	 * Converts the particle deposits into a useful matrix formulation.
	 */
	virtual void initEijMatrix(TMatrixD& eij, TVectorD& truthE);

	/*
	 * Utility method to extract the unique number of detected elements.
	 */
	virtual void populateDetElIndex();

	virtual TVectorD& getProjections(const TMatrixD& eij, TVectorD& proj,
			const TVectorD& truthE) const;

	virtual TMatrixD& getHessian(const TMatrixD& eij, TMatrixD& hess,
			const TVectorD& truthE) const;

	/*
	 * Map to convert detector element to array row/column index.
	 */
	std::map<DetectorElement*, unsigned> myDetElIndex;

};
}

#endif /*LINEARCALIBRATOR_HH_*/
