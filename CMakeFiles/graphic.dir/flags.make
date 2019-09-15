# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# compile CXX with /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++
CXX_FLAGS =   -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk   -std=gnu++11

CXX_DEFINES = -DvtkDomainsChemistry_AUTOINIT="1(vtkDomainsChemistryOpenGL2)" -DvtkIOExport_AUTOINIT="1(vtkIOExportOpenGL2)" -DvtkRenderingContext2D_AUTOINIT="1(vtkRenderingContextOpenGL2)" -DvtkRenderingCore_AUTOINIT="3(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingOpenGL2)" -DvtkRenderingOpenGL2_AUTOINIT="1(vtkRenderingGL2PSOpenGL2)" -DvtkRenderingVolume_AUTOINIT="1(vtkRenderingVolumeOpenGL2)"

CXX_INCLUDES = -I/Users/edward/VTK-8.1.2/Utilities/KWIML -I/Users/edward/VTK-8.1.2/Common/Core -I/Users/edward/VTK-8.1.2/Common/Math -I/Users/edward/VTK-8.1.2/Common/Misc -I/Users/edward/VTK-8.1.2/Common/System -I/Users/edward/VTK-8.1.2/Common/Transforms -I/Users/edward/VTK-8.1.2/Common/DataModel -I/Users/edward/VTK-8.1.2/Common/Color -I/Users/edward/VTK-8.1.2/Common/ExecutionModel -I/Users/edward/VTK-8.1.2/Common/ComputationalGeometry -I/Users/edward/VTK-8.1.2/Filters/Core -I/Users/edward/VTK-8.1.2/Filters/General -I/Users/edward/VTK-8.1.2/Imaging/Core -I/Users/edward/VTK-8.1.2/Imaging/Fourier -I/Users/edward/VTK-8.1.2/ThirdParty/alglib -I/Users/edward/VTK-8.1.2/Filters/Statistics -I/Users/edward/VTK-8.1.2/Filters/Extraction -I/Users/edward/VTK-8.1.2/Infovis/Core -I/Users/edward/VTK-8.1.2/Filters/Geometry -I/Users/edward/VTK-8.1.2/Filters/Sources -I/Users/edward/VTK-8.1.2/Rendering/Core -I/Users/edward/VTK-8.1.2/ThirdParty/zlib -I/Users/edward/VTK-8.1.2/ThirdParty/freetype -I/Users/edward/VTK-8.1.2/Rendering/FreeType -I/Users/edward/VTK-8.1.2/Rendering/Context2D -I/Users/edward/VTK-8.1.2/Charts/Core -I/Users/edward/VTK-8.1.2/Utilities/DICOMParser -I/Users/edward/VTK-8.1.2/ThirdParty/lz4/vtklz4/lib -I/Users/edward/VTK-8.1.2/ThirdParty/lz4/vtklz4 -I/Users/edward/VTK-8.1.2/ThirdParty/lz4 -I/Users/edward/VTK-8.1.2/IO/Core -I/Users/edward/VTK-8.1.2/IO/Legacy -I/Users/edward/VTK-8.1.2/ThirdParty/expat -I/Users/edward/VTK-8.1.2/IO/XMLParser -I/Users/edward/VTK-8.1.2/Domains/Chemistry -I/Users/edward/VTK-8.1.2/Utilities/EncodeString -I/Users/edward/VTK-8.1.2/ThirdParty/glew -I/Users/edward/VTK-8.1.2/Rendering/OpenGL2 -I/Users/edward/VTK-8.1.2/Domains/ChemistryOpenGL2 -I/Users/edward/VTK-8.1.2/IO/XML -I/Users/edward/VTK-8.1.2/Utilities/HashSource -I/Users/edward/VTK-8.1.2/Parallel/Core -I/Users/edward/VTK-8.1.2/Filters/AMR -I/Users/edward/VTK-8.1.2/Filters/FlowPaths -I/Users/edward/VTK-8.1.2/Filters/Generic -I/Users/edward/VTK-8.1.2/Imaging/Sources -I/Users/edward/VTK-8.1.2/Filters/Hybrid -I/Users/edward/VTK-8.1.2/Filters/HyperTree -I/Users/edward/VTK-8.1.2/Imaging/General -I/Users/edward/VTK-8.1.2/Filters/Imaging -I/Users/edward/VTK-8.1.2/Filters/Modeling -I/Users/edward/VTK-8.1.2/Filters/Parallel -I/Users/edward/VTK-8.1.2/Filters/ParallelImaging -I/Users/edward/VTK-8.1.2/Filters/Points -I/Users/edward/VTK-8.1.2/Filters/Programmable -I/Users/edward/VTK-8.1.2/Filters/SMP -I/Users/edward/VTK-8.1.2/Filters/Selection -I/Users/edward/VTK-8.1.2/Filters/Texture -I/Users/edward/VTK-8.1.2/Filters/Topology -I/Users/edward/VTK-8.1.2/ThirdParty/verdict -I/Users/edward/VTK-8.1.2/Filters/Verdict -I/Users/edward/VTK-8.1.2/Utilities/MetaIO/vtkmetaio -I/Users/edward/VTK-8.1.2/Utilities/MetaIO -I/Users/edward/VTK-8.1.2/ThirdParty/jpeg -I/Users/edward/VTK-8.1.2/ThirdParty/png -I/Users/edward/VTK-8.1.2/ThirdParty/tiff/vtktiff/libtiff -I/Users/edward/VTK-8.1.2/ThirdParty/tiff -I/Users/edward/VTK-8.1.2/IO/Image -I/Users/edward/VTK-8.1.2/Imaging/Hybrid -I/Users/edward/VTK-8.1.2/Infovis/Layout -I/Users/edward/VTK-8.1.2/Interaction/Style -I/Users/edward/VTK-8.1.2/Imaging/Color -I/Users/edward/VTK-8.1.2/Rendering/Annotation -I/Users/edward/VTK-8.1.2/Rendering/Volume -I/Users/edward/VTK-8.1.2/Interaction/Widgets -I/Users/edward/VTK-8.1.2/Views/Core -I/Users/edward/VTK-8.1.2/ThirdParty/libproj4/vtklibproj4 -I/Users/edward/VTK-8.1.2/ThirdParty/libproj4 -I/Users/edward/VTK-8.1.2/Geovis/Core -I/Users/edward/VTK-8.1.2/ThirdParty/hdf5/vtkhdf5 -I/Users/edward/VTK-8.1.2/ThirdParty/hdf5 -I/Users/edward/VTK-8.1.2/IO/AMR -I/Users/edward/VTK-8.1.2/IO/EnSight -I/Users/edward/VTK-8.1.2/ThirdParty/netcdf/vtknetcdf/include -I/Users/edward/VTK-8.1.2/ThirdParty/netcdf/vtknetcdf -I/Users/edward/VTK-8.1.2/ThirdParty/netcdf -I/Users/edward/VTK-8.1.2/ThirdParty/exodusII -I/Users/edward/VTK-8.1.2/IO/Exodus -I/Users/edward/VTK-8.1.2/ThirdParty/gl2ps -I/Users/edward/VTK-8.1.2/Rendering/GL2PSOpenGL2 -I/Users/edward/VTK-8.1.2/ThirdParty/libharu/vtklibharu/include -I/Users/edward/VTK-8.1.2/ThirdParty/libharu -I/Users/edward/VTK-8.1.2/IO/Export -I/Users/edward/VTK-8.1.2/IO/ExportOpenGL2 -I/Users/edward/VTK-8.1.2/IO/Geometry -I/Users/edward/VTK-8.1.2/IO/Import -I/Users/edward/VTK-8.1.2/ThirdParty/libxml2/vtklibxml2 -I/Users/edward/VTK-8.1.2/ThirdParty/libxml2 -I/Users/edward/VTK-8.1.2/IO/Infovis -I/Users/edward/VTK-8.1.2/IO/LSDyna -I/Users/edward/VTK-8.1.2/IO/MINC -I/Users/edward/VTK-8.1.2/ThirdParty/oggtheora -I/Users/edward/VTK-8.1.2/IO/Movie -I/Users/edward/VTK-8.1.2/ThirdParty/netcdfcpp -I/Users/edward/VTK-8.1.2/IO/NetCDF -I/Users/edward/VTK-8.1.2/IO/PLY -I/Users/edward/VTK-8.1.2/ThirdParty/jsoncpp -I/Users/edward/VTK-8.1.2/IO/Parallel -I/Users/edward/VTK-8.1.2/IO/ParallelXML -I/Users/edward/VTK-8.1.2/ThirdParty/sqlite -I/Users/edward/VTK-8.1.2/IO/SQL -I/Users/edward/VTK-8.1.2/IO/TecplotTable -I/Users/edward/VTK-8.1.2/IO/Video -I/Users/edward/VTK-8.1.2/Imaging/Math -I/Users/edward/VTK-8.1.2/Imaging/Morphological -I/Users/edward/VTK-8.1.2/Imaging/Statistics -I/Users/edward/VTK-8.1.2/Imaging/Stencil -I/Users/edward/VTK-8.1.2/Interaction/Image -I/Users/edward/VTK-8.1.2/Rendering/ContextOpenGL2 -I/Users/edward/VTK-8.1.2/Rendering/Image -I/Users/edward/VTK-8.1.2/Rendering/LOD -I/Users/edward/VTK-8.1.2/Rendering/Label -I/Users/edward/VTK-8.1.2/Rendering/VolumeOpenGL2 -I/Users/edward/VTK-8.1.2/Views/Context2D -I/Users/edward/VTK-8.1.2/Views/Infovis -isystem /Users/edward/VTK-8.1.2/Utilities/KWSys -isystem /Users/edward/VTK-8.1.2/ThirdParty/hdf5/vtkhdf5/hl/src -isystem /Users/edward/VTK-8.1.2/ThirdParty/hdf5/vtkhdf5/src -isystem /Users/edward/VTK-8.1.2/ThirdParty/netcdfcpp/vtknetcdfcpp 

