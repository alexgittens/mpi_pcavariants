version := "0.0.1"
scalaVersion := "2.10.5"
libraryDependencies ++= Seq(
    "org.scalanlp" %% "breeze" % "0.12",
    "org.scalanlp" %% "breeze-natives" % "0.12"
)

resolvers += "Sonatype Releases" at "https://oss.sonatype.org/content/repositories/releases/"

lazy val submit = taskKey[Unit]("Compute CSFR-0 3D EOFs")
submit <<= (assembly in Compile) map {
  (jarFile: File) => s"src/computeEOFs.sh ${jarFile}" !
} 

val infname = "eofs.hdf5"
val metadatadir = "/global/cscratch1/sd/gittens/conversion-code/ocean_conversion/testOutputs/"
val outfname = "eofs.nc"

lazy val convertSVD = taskKey[Unit]("Convert 2D HDF5 SVD to 3D netCDF")
convertSVD <<= (assembly in Compile) map {
  (jarFile : File) => s"java -jar ${jarFile} $infname $metadatadir $outfname" !
}
