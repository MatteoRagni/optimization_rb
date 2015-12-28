#!/usr/bin/env ruby

# This library is used to generate optimization problem
#
# Author:: Matteo Ragni
# Copyright:: Copyright (2015) University of Trento
# License:: Distributes under the same terms as Ruby

gems = %w|nmatrix nmatrix/lapacke colorize pry|
gems.each do |gem|
  #raise LoadError, "Cannot Load Gem :: #{gem}" unless require gem
  require gem
end

##
# Full container for a set of optimization problem completely defined i Ruby.
# Even if it is not convenient in term of computational time. This is for
# sure a very confortable way to prototype algorigthms.
module Optimization
  load './functions.rb'

  ##
  # This class is a general container for the algorithm that I will define
  class Algorithm
    # Objective function
    attr_reader :objective
    # List of constraints for the problem
    attr_reader :constraints
    # Lagrange functions as a proc
    attr_reader :lagrange, :lagrange_x, :lagrange_xx, :lagrange_lambda
    # History of the optimization problem
    attr_reader :history
    # Options for the algorithm
    attr_accessor :options
    # Size for problem
    attr_reader :size
    # Number of equalities and inequalities
    attr_reader :eq_size, :ineq_size

    ##
    # Initialize a new optimization problem. As for now, multiple dimensions
    # linear functions are not supported
    def initialize(options, objective)
      raise ArgumentError, "options must be an Hash" unless options.is_a? Hash
      raise ArgumentError, "Objective function must be of type ObjectiveFunction" unless
        (objective.is_a? ObjectiveFunction or
         objective.is_a? LinearObjectiveFunction or
         objective.is_a? QuadraticObjectiveFunction)
      @options = options
      @objective = objective
      @size = objective.size
      @constraints = []
      @eq_size = 0
      @ineq_size = 0
    end

    ##
    # Add a new constraint to the problem
    def add_constraint(c)
      raise ArgumentError, "Argument must be a constraint function" unless
        (c.is_a? ConstraintFunction or
         c.is_a? LinearConstraintFunction or
         c.is_a? QuadraticConstraintFunction)
      raise ArgumentError, "Function must accept x of size #{@size}. Is of size #{c.size}" unless @size = c.size
      @eq_size += 1 if c.type == :equality
      @ineq_size += 1 if c.type == :inequality
      @constraints << c
    end

    ##
    # Solves the optimization problem given a starting guess
    def solve(x0)
      raise RuntimeError, "This method is not implemented!"
    end

    def feasible?(x)
      @constraints.each do |c|
        case c.type
        when :equality
          return false if c.f(x) != 0.0
        when :inequality
          return false if c.f(x) < 0.0
        end
      end
    end

    private
    ##
    # Returns the Lagrange function value for a defined problem, in the form of a `Proc`
    def lagrange(x)
      return @lagrange.call(x)
    end

    ##
    # Returns the Lagrange derivative for a defined problem, in the form of a `Proc`
    def lagrange_x(x)
      return @lagrange_x.call(x)
    end

    ##
    # Returns the Lagrange second derivative for a defined problem, in the form of a `Proc`
    def lagrange_xx(x)
      return @lagrange_xx.call(x)
    end

    ##
    # Returns the Lagrange constraints value for a defined problem, in the form of a `Proc`
    def lagrange_lambda(x)
      return @lagrange_lambda.call(x)
    end

    ##
    # This function re-evaluate lagrange function (this will keep lagrange call faster)
    def update_lagrange
      @lagrange = Proc.new do |x|
        value = @objective.f(x)
        @contraints.each do |c|
          value -= c.lamda * c.f(x)
        end
        return value
      end

      @lagrange_x = Proc.new do |x|
        value = @objective.g(x)
        @constraints.each do |c|
          value -= c.lambda * c.g(x)
        end
      end

      @lagrange_xx = Proc.new do |x|
        value = @objective.h(x)
        @constraints.each do |c|
          value -= c.lambda * c.h(x)
        end
      end

      @lagrange_lambda = Proc.new do |x|
        value = NMatrix.zeroes([@constraints.size, 1])
        @constraints.each_with_index do |c, i|
          value[i, 0] = c.f(x)
        end
      end
    end
  end

  ##
  # Class that represents a quadratic optimizer with active set search.
  class QuadraticOptimizer < Algorithm
    # Active set
    attr_reader :active_set

    ##
    # Initializer for our quadratic optimizer with active set
    def initialize(options, objective)
      raise ArgumentError, "Objective function must be a QuadraticObjectiveFunction" unless objective.is_a? QuadraticObjectiveFunction
      super(options, objective)
    end

    ##
    # Add a constraint. It could be only of the linear type
    def add_constraint(c)
      raise ArgumentError, "This algorithm accepts only linear constraint" unless c.is_a? LinearConstraintFunction
      super(c)
    end

    ##
    # Solves algorithm when a
    def solve(x0)
      reorder_constraints
      cholesky
      x = x0
      active_set(x)

      loop do
        z, lambdas = q_problem
        if feasible? z
          x = line_problem(z, x)
          active_set(x)
        else
          x = z
          feasible_set(x)
          @feasible_set.each #XXX

        end

      end
    end

    private
    ##
    # Evaluates the active set for the defined step. Also evaluate active inequality set
    def active_set(x)
      @active_set = []
      @constraints.each_with_index do |c, j|
        case c.type
          when :equality
            @active_set << j
          when :inequality
            @active_set << j if c.active? x
        end
      end
    end

    ##
    # Extracts matrix A and vector b given the active set
    def active_set_matrix
      a_ary = []
      b_ary = []
      @active_set.each do |j|
        a_ary << @constraints[j].a.to_flat_array
        b_ary << @constraints[j].b.to_flat_array
      end
      a_matrix = NMatrix.new([@size, @active_set], a_ary.flatten).transpose
      b_matrix = NMatrix.new([@active_set.size, 1], b_ary.flatten)
      return a_matrix, b_matrix
    end

    ##
    # Evaluate Cholesky decomposition, used in subproblem solution
    def cholesky
      @c, @ct = @objective.s.factorize_cholesky
      @l = @c.invert
      @lt = @l.transpose
      @g = @objective.g
      @d = (@lt.dot(@l)).dot(@g)
      @hinv = @lt.dot(@lt)
    end

    ##
    # Move to the first positions equality constraints. Used to define the active set.
    def reorder_constraints
      eq = []
      ineq = []
      @constraints.each do |c|
        eq << c.dup if c.type == :equality
        ineq << c.dup if c.type == :inequality
      end
      @constraints = []
      eq.each do |c|; @constraints << c; end
      ineq.each do |c|; @constraints << c; end
    end

    ##
    # Returns solution for the sub-problem
    def q_problem
      a, b    = active_set_matrix
      g       = a.transpose.dot(@hinv.dot(a))
      gc, gct = g.factorize_cholesky
      gl      = gc.invert
      glt     = gl.transpose
      ginv    = glt.dot(gl)

      lambdas = ginv.dot(a.transpose.dot(@d) - b)
      xs      = @hinv.dot(a.dot(lambdas) - @g)

      return xs, lambdas
    end
  end
end

if $0 == __FILE__ then
  # Definition
  a = N[[1.0,2.0,3.0],
        [4.0,5.0,6.0],
        [7.0,8.0,9.0]]
  b = N[[-30.0],
        [-20.0],
        [-10.0]]
  c = 35.0
  obj = Optimization::QuadraticObjectiveFunction.new(a, b, c)

  alg = Optimization::QuadraticOptimizer.new({}, obj)
  binding.pry
end
