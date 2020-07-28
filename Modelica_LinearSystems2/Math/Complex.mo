within Modelica_LinearSystems2.Math;
operator record Complex = .Complex "Record defining a Complex number"
  annotation (
    obsolete = "Obsolete operator record - use top-level operator record Complex instead",
    Documentation(info="<html>
<p>
This record extends from the top-level operator record Complex.
It is done for compatibility reasons only and the operator
record will be removed in a&nbsp;future version.
</p>

<p>
This means, generally, that an import such as
</p>
<blockquote><pre>
<strong>import</strong> Modelica_LinearSystems2.Math.Complex;
</pre></blockquote>
<p>
shall be replaced by
</p>
<blockquote><pre>
<strong>import</strong> Complex;
</pre></blockquote>
<p>
in your models.
This new import can even be removed completely in most cases.
</p>
</html>"));
