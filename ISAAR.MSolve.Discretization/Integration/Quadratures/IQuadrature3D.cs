﻿using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Integration.Points;

namespace ISAAR.MSolve.Discretization.Integration.Quadratures
{
	/// <summary>
	/// Collection of integration points that are generated by a traditional 3D quadrature rule, independent of the 
	/// element type. All integration points are with respect to a natural (element local) coordinate system.
	/// These integration points are stored as static fields of an enum class, so that accessing them is fast 
	/// and there is only one copy for all elements. The <see cref="IQuadrature3D"/> object and the collection of integration
	/// points it returns are both immutable.
	/// </summary>
	public interface IQuadrature3D
	{
		/// <summary>
		/// The integration points are sorted based on an order strictly defined for each quadrature.
		/// </summary>
		IReadOnlyList<GaussPoint3D> IntegrationPoints { get; }
	}
}