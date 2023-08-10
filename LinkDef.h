#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ namespace genfit;

// These need no special treatment.
#pragma link C++ class genfit::EventDisplay+;

// these need no special tratment
#pragma link C++ class genfit::AbsKalmanFitter+;
#pragma link C++ class genfit::KalmanFitStatus;
#pragma link C++ class genfit::KalmanFitterRefTrack+;

// these inherit from classes that need custom streamers
#pragma link C++ class genfit::KalmanFittedStateOnPlane+;
#pragma link C++ class genfit::ReferenceStateOnPlane+;

// Classes that needed manually written streamers:
#pragma link C++ class genfit::KalmanFitter-;
#pragma link C++ class genfit::KalmanFitterInfo-;
#pragma link C++ class genfit::DAF-;

#pragma link C++ class genfit::GFGbl+;
#pragma link C++ class genfit::GblFitter+;
#pragma link C++ class genfit::ICalibrationParametersDerivatives+;
#pragma link C++ class genfit::GblFitStatus+;
#pragma link C++ class genfit::GblFitterInfo+;
#pragma link C++ class genfit::GblTrackSegmentController+;
#pragma link C++ class gbl::GblData+;
#pragma link C++ class std::vector<gbl::GblData>+;

#pragma link C++ class genfit::GFRaveVertex+;
#pragma link C++ class genfit::GFRaveTrackParameters+;

#pragma link C++ class genfit::HMatrixU+;
#pragma link C++ class genfit::HMatrixUnit+;
#pragma link C++ class genfit::HMatrixV+;
#pragma link C++ class genfit::HMatrixUV+;
#pragma link C++ class genfit::ProlateSpacepointMeasurement+;
#pragma link C++ class genfit::WireMeasurement+;
#pragma link C++ class genfit::WireMeasurementNew+;
#pragma link C++ class genfit::WirePointMeasurement+;

#pragma link C++ class genfit::HMatrixPhi-;
#pragma link C++ class genfit::FullMeasurement-;
#pragma link C++ class genfit::PlanarMeasurement-;
#pragma link C++ class genfit::SpacepointMeasurement-;

#pragma link C++ class genfit::WireTrackCandHit+;

#pragma link C++ class genfit::RKTrackRep-;
#pragma link C++ class genfit::MplTrackRep-;

#pragma link C++ class genfit::HelixTrackModel+;
#pragma link C++ class genfit::MeasurementCreator+;
#pragma link C++ class genfit::mySpacepointDetectorHit+;
#pragma link C++ class genfit::mySpacepointMeasurement+;

