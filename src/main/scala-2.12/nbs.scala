/**
  * Created by mikesh on 24/11/16.
  */

import com.milaboratory.core.sequence.AminoAcidSequence
import com.milaboratory.core.tree.{SequenceTreeMap, TreeSearchParameters}

import scala.collection.JavaConverters._
import scala.collection.immutable.HashSet
import scala.io.Source

val files = new java.io.File("./vcs/scala-nbs/").listFiles(_.getName.startsWith("top10000"))

val mask0 = Vector.fill(files.length)(false)

val filesWithMask = files.zipWithIndex.map { case (f, i) =>
  (f, mask0.updated(i, true))
}

def mergeMask(m1: Vector[Boolean], m2: Vector[Boolean]): Vector[Boolean] =
  (m1, m2).zipped.map(_ || _)

val stm = new SequenceTreeMap[AminoAcidSequence, (AminoAcidSequence, Int)](AminoAcidSequence.ALPHABET)

val data = filesWithMask.map { case (f, m) =>
  HashSet(Source.fromFile(f).getLines.toList: _*) // unique lines only
    .map { str => new AminoAcidSequence(str) -> m } // convert to milib & assign mask
    .toMap
}.reduceLeft { case (r, m) => // merge datasets
  m.foldLeft(r) { case (dict, (k, v)) =>
    dict + (k -> mergeMask(v, dict.getOrElse(k, mask0)))
  }
}.zipWithIndex // index for converting to matrices

// put to sequence tree
data.foreach { case ((k, _), i) =>
  stm.put(k, (k, i)) // have to put additional copy of sequence to value as they'll be not stored in a "normal" way
}

val tsp = new TreeSearchParameters(1, 0, 0)

val graph = data.flatMap { case ((aa, _), i) =>
  stm.getNeighborhoodIterator(aa, tsp).toList.asScala // best way to query STM
    .filter { case (aa2, _) => !aa.equals(aa2) } // remove exact matches
    .map { case (aa2, j) => if (aa.compareTo(aa2) > 0) (i, j) else (j, i) } // no loops
}.toSet // unique

// Transform graph to sparse matrix

