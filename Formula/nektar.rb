class Nektar < Formula
  desc "High-performance spectral/hp element framework"
  homepage "https://www.nektar.info/"
  url "https://gitlab.nektar.info/nektar/nektar/-/archive/v5.9.0/nektar-v5.9.0.tar.bz2"
  sha256 "8fbb3d96546c72ea5efefe651d9dd921196d075da6c01d5c1c7ddbf964ebd2fe"

  depends_on "arpack"
  depends_on "boost"
  depends_on "cmake"
  depends_on "fftw"
  depends_on "hdf5-mpi"
  depends_on "mdave/nektar/nektar-tinyxml"
  depends_on "numpy"
  depends_on "open-mpi"
  depends_on "opencascade"
  depends_on "pybind11"
  depends_on "python-setuptools"
  depends_on "python@3.14"
  depends_on "scotch"
  depends_on "zlib"

  # Patch for HDF 1.12
  patch :DATA

  def install
    args = std_cmake_args + ["-DNEKTAR_BUILD_TESTS=OFF",
                             "-DNEKTAR_BUILD_UNIT_TESTS=OFF",
                             "-DBoost_NO_WARN_NEW_VERSIONS=1",
                             "-DBoost_USE_MULTITHREADED=OFF",
                             "-DBOOST_ROOT=#{Formula["boost"].opt_prefix}",
                             "-DZLIB_ROOT=#{Formula["zlib"].opt_prefix}"]

    args << "-DNEKTAR_BUILD_DEMOS=OFF"
    args << "-DNEKTAR_BUILD_MESHGEN=ON"
    args << "-DNEKTAR_BUILD_PYTHON=ON"
    args << "-DNEKTAR_SOLVER_DIFFUSION=OFF"
    args << "-DNEKTAR_SOLVER_DUMMY=OFF"
    args << "-DNEKTAR_SOLVER_ELASTICITY=OFF"
    args << "-DNEKTAR_SOLVER_MMF=OFF"
    args << "-DNEKTAR_USE_MPI=ON"
    args << "-DNEKTAR_USE_ARPACK=ON"
    args << "-DNEKTAR_USE_FFTW=ON"
    args << "-DNEKTAR_USE_SCOTCH=ON"
    args << "-DNEKTAR_USE_ARPACK=ON"
    args << "-DNEKTAR_USE_HDF5=ON"
    args << "-DPYTHON_EXECUTABLE=#{Formula["python@3.14"].opt_bin}/python3.14"

    mkdir "build" do
      system "cmake", "..", *args
      system "make", "install"

      # Also need to install NekPy bindings
      python = Formula["python@3.14"].opt_bin/"python3.14"
      cd "python" do
        system python, "-m", "pip", "install", *std_pip_args(prefix: libexec), "."
      end

      site_packages = Language::Python.site_packages(python)
      pth_contents = "import site; site.addsitedir('#{libexec/site_packages}')\n"
      (prefix/site_packages/"homebrew-nektar.pth").write pth_contents

      Dir[libexec/site_packages/"NekPy/**/*.so"].each do |so|
        MachO::Tools.add_rpath(so, lib.to_s)
      end
    end
  end

  test do
    (testpath/"helm.cpp").write <<~EOS
      #include <iostream>
      #include <MultiRegions/ContField.h>
      #include <SpatialDomains/MeshGraph.h>
      #include <SpatialDomains/MeshGraphIO.h>

      using namespace Nektar;
      using namespace std;

      int main(int argc, char *argv[])
      {
          LibUtilities::SessionReaderSharedPtr vSession
              = LibUtilities::SessionReader::CreateInstance(argc, argv);
          SpatialDomains::MeshGraphSharedPtr graph =
              SpatialDomains::MeshGraphIO::Read(vSession);
          MultiRegions::ContFieldSharedPtr fld =
              MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
                  vSession, graph, vSession->GetVariable(0));

          // Set up forcing
          const int nq = fld->GetNpoints(), nm = fld->GetNcoeffs();
          Array<OneD, NekDouble> xc0(nq), xc1(nq), xc2(nq);
          fld->GetCoords(xc0, xc1, xc2);
          vSession->GetFunction("Forcing", 0)->Evaluate(xc0, xc1, xc2, fld->UpdatePhys());

          // Solve Helmholtz
          StdRegions::ConstFactorMap factors;
          factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
          Vmath::Zero(nm, fld->UpdateCoeffs(), 1);
          fld->HelmSolve(fld->GetPhys(), fld->UpdateCoeffs(), factors);

          // Evaluate L^2 error
          fld->BwdTrans(fld->GetCoeffs(), fld->UpdatePhys());
          Array<OneD, NekDouble> sol(nq);
          vSession->GetFunction("ExactSolution", 0)->Evaluate(xc0, xc1, xc2, sol);
          cout << fld->L2(fld->GetPhys(), sol) << endl;
          return 0;
      }
    EOS

    (testpath/"CMakeLists.txt").write <<~EOS
      cmake_minimum_required(VERSION 4.0)
      set(CMAKE_CXX_STANDARD 17)
      find_package(Nektar++ REQUIRED NO_MODULE NO_DEFAULT_PATH NO_CMAKE_BUILDS_PATH NO_CMAKE_PACKAGE_REGISTRY)
      find_package(HDF5)
      include_directories(${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
      add_definitions(${NEKTAR++_DEFINITIONS})
      link_directories(${NEKTAR++_LIBRARY_DIRS} ${NEKTAR++_TP_LIBRARY_DIRS})
      add_executable(helm helm.cpp)
      target_link_libraries(helm ${NEKTAR++_LIBRARIES})
    EOS

    (testpath/"project.py").write <<~EOS
      from NekPy.LibUtilities import PointsKey, PointsType, BasisKey, BasisType
      from NekPy.StdRegions import StdQuadExp
      import numpy as np
      ptsKey   = PointsKey(9, PointsType.GaussLobattoLegendre)
      basisKey = BasisKey(BasisType.Modified_A, 8, ptsKey)
      quadExp  = StdQuadExp(basisKey, basisKey)
      x, y     = quadExp.GetCoords()
      fx       = np.sin(x) * np.cos(y)
      proj     = quadExp.FwdTrans(fx)
      assert np.allclose(fx, quadExp.BwdTrans(proj))
    EOS

    (testpath/"input.xml").write <<~EOS
      <NEKTAR>
          <GEOMETRY DIM="2" SPACE="2">
              <VERTEX COMPRESSED="B64Z-LittleEndian" BITSIZE="64">eJxjYEAGD/aj0gwMjAzYAEKeCas8AjDj0AcDLNjNt4exWLG7Dy7PhlUfwh52HObCAAd2/XB1AAerEJkA</VERTEX>
              <EDGE COMPRESSED="B64Z-LittleEndian" BITSIZE="64">eJxjYEAFjDhoJhw0Mw4aBljQ1MP4rDj4bGh8mHnsaO6BqeNA48PUcaLxYfZzoYnD9HOj8WHuAgBAqACe</EDGE>
              <ELEMENT>
                  <Q COMPRESSED="B64Z-LittleEndian" BITSIZE="64">eJxjYEAFjFCaCUoz4xBngdKsUJoNTZ4dSnNAaU40c5jRxLmgNDea+QAUwABZ</Q>
              </ELEMENT>
              <COMPOSITE>
                  <C ID="0"> Q[0-3] </C>
                  <C ID="1"> E[0,3,5,6,7,8,10,11] </C>
              </COMPOSITE>
              <DOMAIN> C[0] </DOMAIN>
          </GEOMETRY>
          <EXPANSIONS>
              <E COMPOSITE="C[0]" NUMMODES="9" TYPE="MODIFIED" FIELDS="u" />
          </EXPANSIONS>
          <CONDITIONS>
              <PARAMETERS> <P> Lambda = 1 </P> </PARAMETERS>
              <VARIABLES> <V ID="0"> u </V> </VARIABLES>
              <BOUNDARYREGIONS> <B ID="0"> C[1] </B> </BOUNDARYREGIONS>
              <BOUNDARYCONDITIONS>
                  <REGION REF="0">
                      <D VAR="u" VALUE="sin(PI*x)*sin(PI*y)" />
                  </REGION>
              </BOUNDARYCONDITIONS>
              <FUNCTION NAME="Forcing">
                  <E VAR="u" VALUE="-(Lambda + 2*PI*PI)*sin(PI*x)*sin(PI*y)" />
              </FUNCTION>
              <FUNCTION NAME="ExactSolution">
                  <E VAR="u" VALUE="sin(PI*x)*sin(PI*y)" />
              </FUNCTION>
          </CONDITIONS>
      </NEKTAR>
    EOS
    system "cmake", "-DNektar++_DIR=#{lib}/nektar++/cmake/", "."
    system "make", "VERBOSE=1"
    assert (`./helm input.xml`.to_f < 1.0e-8)
    python = Formula["python@3.14"]
    system python.bin/"python3.14", "project.py"
  end
end

__END__
diff --git a/cmake/ThirdPartyPython.cmake b/cmake/ThirdPartyPython.cmake
index f1e37ba1a..643bf225c 100644
--- a/cmake/ThirdPartyPython.cmake
+++ b/cmake/ThirdPartyPython.cmake
@@ -40,18 +40,27 @@ IF (NEKTAR_BUILD_PYTHON)
         COMMAND ${Python3_EXECUTABLE} -m pip install .
         WORKING_DIRECTORY ${NEKPY_BASE_DIR})
 
-    EXTERNALPROJECT_ADD(
-        pybind11
-        PREFIX ${TPSRC}
-        URL ${TPURL}/pybind11-002c05b17.zip
-        URL_MD5 4fd2f2df6a8a4f28c260e56163bf4481
-        STAMP_DIR ${TPBUILD}/stamp
-        DOWNLOAD_DIR ${TPSRC}
-        SOURCE_DIR ${TPSRC}/pybind11
-        BINARY_DIR ${TPBUILD}/pybind11
-        TMP_DIR ${TPBUILD}/pybind11-tmp
-        INSTALL_DIR ${TPDIST}
-        CONFIGURE_COMMAND ${CMAKE_COMMAND}
+    # Require pybind11 v3 or later for unique_ptr transfer
+    FIND_PACKAGE(pybind11 3 CONFIG QUIET)
+
+    IF (pybind11_FOUND)
+        MESSAGE(STATUS "Found pybind11: ${pybind11_INCLUDE_DIR}")
+        ADD_CUSTOM_TARGET(pybind11 ALL)
+        INCLUDE_DIRECTORIES(${pybind11_INCLUDE_DIR})
+    ELSE()
+        MESSAGE(STATUS "Build pybind11: ${TPDIST}/include/pybind11")
+        EXTERNALPROJECT_ADD(
+            pybind11
+            PREFIX ${TPSRC}
+            URL ${TPURL}/pybind11-3.0.1.zip
+            URL_MD5 53f015a45ffaeeec2ad605ac436526ba
+            STAMP_DIR ${TPBUILD}/stamp
+            DOWNLOAD_DIR ${TPSRC}
+            SOURCE_DIR ${TPSRC}/pybind11
+            BINARY_DIR ${TPBUILD}/pybind11
+            TMP_DIR ${TPBUILD}/pybind11-tmp
+            INSTALL_DIR ${TPDIST}
+            CONFIGURE_COMMAND ${CMAKE_COMMAND}
             -G ${CMAKE_GENERATOR}
             -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
             -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
@@ -59,12 +68,13 @@ IF (NEKTAR_BUILD_PYTHON)
             -DPYBIND11_TEST=OFF
             ${TPSRC}/pybind11)
 
-     # Add third-party include to include path.
-     INCLUDE_DIRECTORIES(${TPDIST}/include)
+        # Add third-party include to include path.
+        INCLUDE_DIRECTORIES(${TPDIST}/include)
 
-     ADD_DEPENDENCIES(thirdparty pybind11)
+        ADD_DEPENDENCIES(thirdparty pybind11)
+    ENDIF()
 
-     FILE(WRITE ${NEKPY_BASE_DIR}/NekPy/__init__.py "# Adjust dlopen flags to avoid OpenMPI issues
+    FILE(WRITE ${NEKPY_BASE_DIR}/NekPy/__init__.py "# Adjust dlopen flags to avoid OpenMPI issues
 try:
     import DLFCN as dl
     import sys
diff --git a/cmake/ThirdPartyBoost.cmake b/cmake/ThirdPartyBoost.cmake
index c71c368960..22c962efb7 100644
--- a/cmake/ThirdPartyBoost.cmake
+++ b/cmake/ThirdPartyBoost.cmake
@@ -12,29 +12,29 @@ MESSAGE(STATUS "Searching for Boost:")
 # Minimum version and boost libraries required
 SET(MIN_VER "1.60.0")
 IF (NEKTAR_USE_BOOST_FILESYSTEM)
-    SET(NEEDED_BOOST_LIBS iostreams filesystem system program_options)
+    SET(NEEDED_BOOST_LIBS iostreams filesystem program_options)
 ELSE()
-    SET(NEEDED_BOOST_LIBS iostreams system program_options)
+    SET(NEEDED_BOOST_LIBS iostreams program_options)
 ENDIF()
 
 SET(Boost_NO_BOOST_CMAKE ON)
 IF( BOOST_ROOT )
     SET(Boost_NO_SYSTEM_PATHS ON)
-    FIND_PACKAGE( Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS})
+    FIND_PACKAGE( Boost ${MIN_VER} COMPONENTS ${NEEDED_BOOST_LIBS} OPTIONAL_COMPONENTS system)
 ELSE ()
     SET(TEST_ENV1 $ENV{BOOST_HOME})
     SET(TEST_ENV2 $ENV{BOOST_DIR})
     IF (DEFINED TEST_ENV1)
         SET(BOOST_ROOT $ENV{BOOST_HOME})
         SET(Boost_NO_SYSTEM_PATHS ON)
-        FIND_PACKAGE( Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS} )
+        FIND_PACKAGE( Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS} OPTIONAL_COMPONENTS system)
     ELSEIF (DEFINED TEST_ENV2)
         SET(BOOST_ROOT $ENV{BOOST_DIR})
         SET(Boost_NO_SYSTEM_PATHS ON)
-        FIND_PACKAGE( Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS} )
+        FIND_PACKAGE( Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS} OPTIONAL_COMPONENTS system)
     ELSE ()
         SET(BOOST_ROOT ${TPDIST})
-        FIND_PACKAGE( Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS} )
+        FIND_PACKAGE( Boost ${MIN_VER} QUIET COMPONENTS ${NEEDED_BOOST_LIBS} OPTIONAL_COMPONENTS system)
     ENDIF()
 ENDIF()
 

diff --git a/cmake/ThirdPartyMPI.cmake b/cmake/ThirdPartyMPI.cmake
index 2f538c297e..8c04b741fb 100644
--- a/cmake/ThirdPartyMPI.cmake
+++ b/cmake/ThirdPartyMPI.cmake
@@ -77,6 +77,7 @@ IF( NEKTAR_USE_MPI )
                 -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                 -DCMAKE_BUILD_TYPE:STRING=Debug
                 -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
+                -DCMAKE_POLICY_VERSION_MINIMUM=3.5
                 ${TPSRC}/gsmpi-1.2.1_2
         )
         THIRDPARTY_LIBRARY(GSMPI_LIBRARY STATIC gsmpi DESCRIPTION "GSMPI Library")

