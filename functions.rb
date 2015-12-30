#!/usr/bin/env ruby
require 'pry'
# This library contains a series of helper functions, for example
#
#  * Non linear function defined from `n` to `1`
#  * Linear function in form `A * x + b` defined from `n` to `m`
#  * Quadratic function defined from `n` to `1`
#  * Constraints, as `:equality` and `:inequality`
#
# Author:: Matteo Ragni
# Copyright:: Copyright (2015) University of Trento
# License:: Distributes under the same terms as Ruby

##
# Loads NMatrix if not already loaded
begin
  NMatrix
rescue NameError
  require 'nmatrix'
end

module Optimization
  ##
  # Contains a definition for a general purpose function, that may be used as
  # objective function or as constraint
  class Function
    # Input size (input check not enabled due to performance)
    attr_reader :size

    ##
    # Initializer for a new function. Takes two input, the size of input of the function
    # and the level of derivative definition
    def initialize(x_size = 1, derivative = :gradient)
      raise ArgumentError, "Size of vector x must be a positive integer" unless x_size.is_a? Fixnum
      raise ArgumentError, "Size of vector x must be a positive integer" unless x_size >= 0
      raise ArgumentError, "Definition of derivative level must be a symbol" unless derivative.is_a? Symbol
      raise ArgumentError, "Definition of derivative level must be :function, :gradient or :hessian" unless [:function, :gradient, :hessian].include? derivative

      @size = x_size
      @f = nil

      if derivative == :gradient || derivative == :hessian
        @g = nil
      else
        @g = Proc.new { |x| NMatrix.zeroes [@size] }
      end

      if derivative == :hessian
        @h = nil
      else
        @h = Proc.new { |x| NMatrix.zeroes [@size, @size] }
      end
    end

    ##
    # Define the function through a Proc
    def f=(p)
      raise ArgumentError, "Argument must be a Proc" unless p.is_a? Proc
      @f = p.dup
    end

    ##
    # Define a gradient through a Proc
    def g=(p)
      raise ArgumentError, "Argument must be a Proc" unless p.is_a? Proc
      @g = p.dup
    end

    ##
    # Define an Hessian through a Proc
    def h=(p)
      raise ArgumentError, "Argument must be a Proc" unless p.is_a? Proc
      @h = p.dup
    end

    ##
    # Call the function. **No trasformation to Float is performed please pay attention**
    def f(x)
      @f.call(x)
    end

    ##
    # Call the gradient. **No trasformation to Float is performed please pay attention**
    def g(x)
      @g.call(x)
    end

    ##
    # Call the Hessian **No trasformation to Float is performed please pay attention**
    def h(x)
      @h.call(x)
    end

    ##
    # Test if vector in input could be used and return results to check if functions
    # are correct. For testing only
    def test(x)
      check_input(x)
      raise ArgumentError, "function not defined" unless @f
      raise ArgumentError, "gradient not defined" unless @g
      raise ArgumentError, "hessian not defined" unless @h

      f_test = f(x)
      g_test = g(x)
      h_test = h(x)
      check_output(f_test, g_test, h_test)

      return { f: f_test, g: g_test, h: h_test }
    end

    private
    ##
    # Check if a x input is correct for our problem
    def check_input(x)
      raise ArgumentError, "x must be a Vector #{@size}x1" unless x.is_a? NMatrix
      raise ArgumentError, "x must be a Vector #{@size}x1. Too many columns" unless x.dim == 1
      raise ArgumentError, "x must be a Vector #{@size}x1. Too many element in vector" unless x.shape[1] == @size
    end

    ##
    # Check if our problem gives a correct output (due to the fact that is a Proc)
    def check_output(f_test, g_test, h_test)
      raise RuntimeError, "function does not return correct type" unless f_test.is_a? Float
      raise RuntimeError, "gradient does not return correct type" unless f_test.is_a? Float
      raise RuntimeError, "function does not return correct type" unless f_test.is_a? Float
    end
  end

  ##
  # Generalized non linear function for an Objective function. Inherits completely
  # from Optimization::Function class
  class ObjectiveFunction < Function; end

  ##
  # Class used to rapresent general non linear constraints in form of
  #  * `:equality` -> f(x) = 0
  #  * `:inequality` -> f(x) >= 0
  class ConstraintFunction < Function
    # Lambda multiplier as defined in optimization litterature
    attr_accessor :lambda
    # A constraint should be of `:equality` or `:inequality` type
    attr_reader :type

    ##
    # Initializer for non linear constraints. Please note that a symbol of type
    # `:equality` or `:inequality`. Other parameters equal to Optimization::Function::new
    def initialize(type, x_size = 1, derivative = :gradient)
      raise ArgumentError, "Type must be a Symbol: :equality or :inequality" unless type.is_a? Symbol
      raise ArgumentError, "Type must be a Symbol: :equality or :inequality" unless [:equality, :inequality].include? type
      super(x_size, derivative)
      @type = type
      @lambda = 0.0
    end

    ##
    # Returns if a constraint is active or not. Please notice
    # that an equality constrain should always be active.
    def active?(x); f(x) <= 0; end
  end

  ##
  # Class used to rapresent a single linear function, in the form `a' * x + b`,
  # in which `a in R(n x 1)` and `b` is a scalar value
  class LinearFunction < Function
    # the linear coefficient vector, should be `n x 1`. Also transposed version is saved
    attr_reader :a, :at
    # the linear intercept scalar
    attr_reader :b
    # Input size is retireved from input vector dimension
    attr_reader :size

    ##
    # Initializer for single linear function. Vector a and b must be provided
    def initialize(a, b)
      raise ArgumentError, "The first input should be a matrix A (A*x + b)" unless a.is_a? NMatrix
      raise ArgumentError, "The second input should be a matrix B (A*x + b)" unless b.is_a? Float
      raise ArgumentError, "a must be of dimensions #{a.shape[0]} x 1" unless a.shape[1] == 1

      @a = a.dup
      @at = @a.transpose.dup
      @b = b
      @size = @a.shape[0]
    end

    ##
    # Return function value
    def f(x)
      (@at.dot(x) + @b)[0]
    end
    ##
    # Return analitycal gradient value
    def g(x)
      @at
    end
    ##
    # Return analitycal hessina value
    def h(x)
      return NMatrix.zeroes([@size, @size])
    end

    def test(x)
      check_input(x)
      return { f: f(x), g: g(x), h: h(x) }
    end

    # Removing functions that will be not used
    [:f=, :g=, :h=, :test_output].each do |m|
      undef_method m if self.instance_methods.include?(m)
    end
  end

  ##
  # Generalized non linear function for an Objective function. Inherits completely
  # from Optimization::Function class
  class LinearObjectiveFunction < LinearFunction; end

  ##
  # Class used to rapresent general linear constraints in form of
  #  * `:equality` -> a' * x + b = 0
  #  * `:inequality` -> a' * x + b >= 0
  class LinearConstraintFunction < LinearFunction
    # Lambda multiplier as defined in optimization litterature
    attr_accessor :lambda
    # A constraint should be of `:equality` or `:inequality` type
    attr_reader :type

    ##
    # Initializer. Refer to Optimization::ConstraintFunction::new
    # and Optimization::LinearFunction::new
    def initialize(type, a, b)
      raise ArgumentError, "Type must be a Symbol: :equality or :inequality" unless type.is_a? Symbol
      raise ArgumentError, "Type must be a Symbol: :equality or :inequality" unless [:equality, :inequality].include? type
      super(a, b)
      @type = type
      @lambda = 0.0
    end

    ##
    # Returns if a constraint is active or not. Please notice
    # that an equality constrain should always be active.
    def active?(x); f(x) <= 0; end

    ##
    # Returns the lagrange term for the constraint
    def lagrange_term(x)
      return @lambda * f(x)
    end
  end


  ##
  # Class used to represent a linear constraint, of multiple dimension, in the form
  # `a x + b` in which `a in R(m x n)` and `b in R(m x 1)`. Due to efficiency
  # reasons, this class is implemented as a series of function.
  class LinearFunctions < Function
    # The linear coefficient matrix, should be `m x n`
    attr_reader :a
    # The linear intercept vector, should be `m x 1`
    attr_reader :b
    # Input size is retrieved from matrix dimensions
    attr_reader :size
    # Output size is retrieved from matrix dimensions
    attr_reader :cosize

    ##
    # Initializer for a linear function. Matrix A and vector B must be provided
    def initialize(a, b)
      raise ArgumentError, "The first input should be a matrix A (A*x + b)" unless a.is_a? NMatrix
      raise ArgumentError, "The second input should be a matrix B (A*x + b)" unless b.is_a? NMatrix
      raise ArgumentError, "a.size[0] must be equal to b.size[0]" unless a.shape[0] == b.shape[0]
      raise ArgumentError, "b must be of dimensions #{b.shape[0]} x 1" unless b.shape[1] == 1

      @a = a.dup
      @b = b.dup
      @size = @a.shape[1]
      @cosize = @a.shape[0]
    end

    ##
    # Return function value
    def f(x)
      @a.dot(x) + @b
    end
    ##
    # Return analitycal gradient value
    def g(x)
      @a
    end
    ##
    # Return analitycal hessina value
    def h(x)
      return NMatrix.zeroes([@cosize, @size, @size])
    end

    def test(x)
      check_input(x)
      return { f: f(x), g: g(x), h: h(x) }
    end

    # Removing functions that will be not used
    [:f=, :g=, :h=, :test_output].each do |m|
      undef_method m if self.instance_methods.include?(m)
    end
  end

  ##
  # Generalized non linear function for an Objective function. Inherits completely
  # from Optimization::Function class
  class LinearObjectiveFunctions < LinearFunctions; end

  ##
  # Class used to rapresent general linear constraints in form of
  #  * `:equality` -> A * x + b = 0
  #  * `:inequality` -> A * x + b >= 0
  class LinearConstraintFunctions < LinearFunctions
    # Lambda multiplier as defined in optimization litterature
    attr_accessor :lambda
    # A constraint should be of `:equality` or `:inequality` type
    attr_reader :type

    ##
    # Initializer. Refer to Optimization::ConstraintFunction::new
    # and Optimization::LinearFunction::new
    def initialize(type, a, b)
      raise ArgumentError, "Type must be a Symbol: :equality or :inequality" unless type.is_a? Symbol
      raise ArgumentError, "Type must be a Symbol: :equality or :inequality" unless [:equality, :inequality].include? type
      super(a, b)
      @type = type
      @lambda = NMatrix.new([1, @cosize], 0.0)
    end

    ##
    # Returns if a constraint is active or not. Please notice
    # that an equality constrain should always be active.
    def active?(x); f(x) <= 0; end

    ##
    # Returns the lagrange term for the constraint
    def lagrange_term(x)
      return @lambda.dot f(x)
    end
  end

  ##
  # This class represents a general quadratic function in the form
  # `1/2 (x' * A * x) + b' * x + c`. Remember that this kind of function returns
  # an output of dimension 1. Consider that:
  #  * `A in R(n x n)`
  #  * `b in R(n x 1)`
  #  * `c in R(1 x 1)`
  # and, once provided the matrix A, this will be decomposed in a S + W, that
  # symmetric and antisymmetric parts: `A = S + W = (A + A')/2 + (A - A')/2`
  class QuadraticFunction < Function
    # The quadratic symmetric coefficient matrix, should be `n x n`
    attr_reader :s
    # The quadratic anti-symmetric coefficient matrix, should be `n x n`
    attr_reader :w
    # The complete matrix A, should be `n x n`
    attr_reader :a
    # The linear coefficient vector, should be `n x 1`. It is automatically transposed
    attr_reader :bt, :b
    # Input size is retrieved from matrix dimensions
    attr_reader :size

    ##
    # Initializer for a linear function. Matrix A and vector B must be provided
    def initialize(a, b, c)
      raise ArgumentError, "The first input should be a matrix A (1/2 (x' * A * x) + b' * x + c)" unless a.is_a? NMatrix
      raise ArgumentError, "The second input should be a matrix b (1/2 (x' * A * x) + b' * x + c)" unless b.is_a? NMatrix
      raise ArgumentError, "The third input should be a Float c (1/2 (x' * A * x) + b' * x + c)" unless c.is_a? Float
      raise ArgumentError, "A.shape[0] must be equal to b.shape[0]" unless a.shape[0] == b.shape[0]
      raise ArgumentError, "A must be square matrix" unless a.shape[0] == a.shape[1]

      @a = a.dup
      @s = (@a + @a.transpose) * 0.5
      @w = (@a - @s)
      @b = b
      @bt = @b.transpose.dup
      @c = N[[c]]
      @size = @a.shape[0]
    end

    ##
    # Returns intercept value, should be `1 x 1`, and it is a method because
    # it is converted automatically as a Float (stored as a NMatrix)
    def c; c.to_f; end

    ##
    # Return function value. In function `(f)[0]` is used to get directly the
    # the value of the function instead of the NMatrix object
    def f(x)
      (x.transpose.dot((@s.dot x)) + @bt.dot(x) + @c)[0]
    end
    ##
    # Return analitycal gradient value
    def g(x)
      x.transpose.dot(@s) + @bt
    end
    ##
    # Return analitycal hessian value
    def h(x)
      @s
    end

    def test(x)
      check_input(x)
      return { f: f(x), g: g(x), h: h(x) }
    end

    # Removing functions that will be not used
    [:f=, :g=, :h=, :test_output].each do |m|
      undef_method m if self.instance_methods.include?(m)
    end
  end

  ##
  # Generalized quadratic function for an Objective function. Inherits completely
  # from Optimization::QuadraticFunction class
  class QuadraticObjectiveFunction < QuadraticFunction; end

  ##
  # Class used to rapresent general quadratic constraints in form of
  #  * `:equality` -> 1/2 x' * A * x + b' * x + c = 0
  #  * `:inequality` -> 1/2 x' * A * x + b' * x + c = 0
  class QuadraticConstraintFunction < QuadraticFunction
    # Lambda multiplier as defined in optimization litterature
    attr_accessor :lambda
    # A constraint should be of `:equality` or `:inequality` type
    attr_reader :type

    ##
    # Initializer. Refer to Optimization::ConstraintFunction::new
    # and Optimization::LinearFunction::new
    def initialize(type, a, b, c)
      raise ArgumentError, "Type must be a Symbol: :equality or :inequality" unless type.is_a? Symbol
      raise ArgumentError, "Type must be a Symbol: :equality or :inequality" unless [:equality, :inequality].include? type
      super(a, b, c)
      @type = type
      @lambda = 0.0
    end

    ##
    # Returns if a constraint is active or not. Please notice
    # that an equality constrain should always be active.
    def active?(x); f(x) <= 0; end

    ##
    # Returns the lagrange term for the constraint
    def lagrange_term(x)
      return @lambda * f(x)
    end
  end
end

## End of Library #########################################################################

##
# LIBRARY TESTING!!
if __FILE__ == $0 then # :nodoc
  require 'pp'

  x = N[[1.0],[2.0],[3.0]] # Vettore colonna 3 x 1

  # As for now the test will be only constraints
  ## NON LINEAR ######################################################
  c_nl = Optimization::ConstraintFunction.new(:inequality, 3, :hessian)
  c_nl.f = Proc.new { |x| x[0]**3 + Math::sin(x[1]) + Math::cos(x[2]) }
  c_nl.g = Proc.new { |x| N[[3.0 * x[0]**2], [-Math::cos(x[1])], [Math::sin(x[2])]] }
  c_nl.h = Proc.new { |x| N[[6*x[0], 0, 0], [0, -Math::sin(x[1]), 0], [0, 0, -Math::cos(x[2])]] }

  puts <<-REPORT
%%%%%%
%% NON LINEAR
Testing nonlinear function #{c_nl.type} constraint:

 in point #{x}

 - size: #{c_nl.size}
 - lambda: #{c_nl.lambda}
 - is active in x? #{c_nl.active?(x)}

REPORT
  print "f(x) = "
  pp c_nl.f(x)
  puts
  print "g(x) = "
  pp c_nl.g(x)
  puts
  print "h(x) = "
  pp c_nl.h(x)
  puts

  ## LINEAR ##########################################################
  a = (N[[1.0,2.0,3.0]]).transpose
  b = 10.0
  c_l = Optimization::LinearConstraintFunction.new(:inequality, a, b)

  puts <<-REPORT
%%%%%%
%% LINEAR
Testing linear function #{c_l.type} constraint:

 in point #{x}

 - size: #{c_l.size}
 - lambda: #{c_l.lambda}
 - is active in x? #{c_l.active?(x)}

REPORT
  puts "A = "
  pp c_l.a
  puts "b = "
  pp c_l.b
  puts

  print "f(x) = "
  pp c_l.f(x)
  puts
  print "g(x) = "
  pp c_l.g(x)
  puts
  print "h(x) = "
  pp c_l.h(x)
  puts


  ## LINEAR ##########################################################
  a = N[[1.0,2.0,3.0],
        [4.0,5.0,6.0],
        [7.0,8.0,9.0]]
  b = N[[-30.0],
        [-20.0],
        [-10.0]]
  c_l = Optimization::LinearConstraintFunctions.new(:inequality, a, b)

  puts <<-REPORT
%%%%%%
%% LINEARS
Testing linear function #{c_l.type} constraint:

 in point #{x}

 - size: #{c_l.size}
 - output size: #{c_l.cosize}
 - lambda: #{c_l.lambda}
 - is active in x? #{c_l.active?(x)}

REPORT
puts "A = "
pp c_l.a
puts "b = "
pp c_l.b
puts

  print "f(x) = "
  pp c_l.f(x)
  puts
  print "g(x) = "
  pp c_l.g(x)
  puts
  print "h(x) = "
  pp c_l.h(x)
  print "lagrange_term(x) = "
  pp c_l.lagrange_term(x)
  puts

  ## QUADRATIC ##########################################################
  a = N[[1.0,2.0,3.0],
        [4.0,5.0,6.0],
        [7.0,8.0,9.0]]
  b = N[[-30.0],
        [-20.0],
        [-10.0]]
  c = 35.0
  c_q = Optimization::QuadraticConstraintFunction.new(:inequality, a, b, c)

  puts <<-REPORT
%%%%%%
%% QUADRATIC
Testing quadratic function #{c_q.type} constraint:

 in point #{x}

 - size: #{c_q.size}
 - lambda: #{c_q.lambda}
 - is active in x? #{c_q.active?(x)}

REPORT
  print "f(x) = "
  pp c_q.f(x)
  puts
  print "g(x) = "
  pp c_q.g(x)
  puts
  print "h(x) = "
  pp c_q.h(x)
  puts

end
