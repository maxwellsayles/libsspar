/**
 * Convert a C file of just primorial representations into
 * one that contains documentations, the representations, and a lookup
 * table for each representation.
 */

val header = """
#include <stdint.h>

/**
 * Representation is given as an array 'A' such that
 * N = ((((A[0].b*A[0].a) \\pm A[1].b)*A[1].a) \\pm A[2].b)*A[2].a ....
 *
 * where A[i].b is added if the high bit of 'b' is clear,
 * and subtracted if the high bit of 'b' is set.
 *
 * We use the high bit so that if 'b' is 0, we can still distinguish between
 * a 'b' that should be added and a 'b' that should be subtracted.
 */
typedef struct {
    uint16_t a;
    uint16_t b;
} factored_two_three_term16_t;

"""

if (args.length < 2) { 
  println("Usage: scala convert.scala <type> <infile> <outfile>")
  sys.exit(0)
}

val qform_type = args(0)

def convert(inFilename: String, outFilename: String) {
  val out = new java.io.FileWriter(outFilename)
  out.write(header)

  val sizeRegex = """primorial(\d+)\[(\d+)\]""".r
  val lines = scala.io.Source.fromFile(inFilename).getLines.toList
  lines foreach { line => out.write(line + "\n") }
  val sizes = lines.foldLeft(Map[Int, Int]()) { (map, line) =>
    sizeRegex findFirstIn line match {
      case Some(sizeRegex(i, size)) =>
	val i2 = i.toInt
	val size2 = size.toInt
	map + (i2 -> size2)
      case None => map
    }
  }

  val minSize = sizes.keys.min
  val maxSize = sizes.keys.max
  val termCounts =
    "0, " * minSize + ((minSize to maxSize) map { sizes(_).toString } mkString ", ")
  val labels =
    "0, " * minSize + ((minSize to maxSize) map { "primorial%d" format _ } mkString ", ")
  val labelsMsg =
    "const factored_two_three_term16_t* %s_primorial_terms[] = {%s};".format(qform_type, labels.tail.tail.toString)
  val termCountsMsg =
    "const uint16_t %s_primorial_term_counts[] = {%s};".format(qform_type, termCounts.tail.tail.toString)

  out.write(labelsMsg + "\n")
  out.write(termCountsMsg + "\n")
  out.close()
}

convert(args(1), args(2))

