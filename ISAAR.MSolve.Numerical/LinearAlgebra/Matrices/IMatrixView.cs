﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra.Commons;
using ISAAR.MSolve.Numerical.LinearAlgebra.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Matrices
{
    public interface IMatrixView: IEntrywiseOperable, IIndexable2D, IWriteable
    {
        IMatrixView MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false);
        IMatrixView MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false);
        IVectorView MultiplyRight(IVectorView vector, bool transposeThis = false);

        IMatrixView Transpose();
    }
}
