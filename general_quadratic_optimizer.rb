module Optimization
  ##
  # Class to perform a simple constrained optimization
  class QuadraticOptimizer < Algorithm
    ##
    # Initializer for quadratic optimization algorithm
    def initialize(options, objective)
      raise ArgumentError, "Objective function must be a QuadraticObjectiveFunction" unless objective.is_a? QuadraticObjectiveFunction
      super options, objective

      # Initially allocated elements
      @h     = @objective.s
      @h_inv = @h.invert
      @g     = @objective.g
      @d     = @hinv.dot @g
    end

    ##
    # Add a constraint to the problem. They could be only equality contraints, so inequalities will be converted in equalities
    def add_constraints(c)
      raise ArgumentError, "This algorithm accepts only linear constraint" unless c.is_a? LinearConstraintFunction
      super c.to_equality
    end

    def solve
      constraint_matrix                                                 # A, b
      @w = (@a.transpose).dot(@hinv.dot(@a))                                  # W = (A^T) H^(-1) A
      @w_pinv = (((@w.traspose).dot(@w)).inverse).dot(@w.transpose)       # W+ = (W^T W)^(-1) W^T

      @lambdas = @w_pinv.dot((@a.traspose).dot(@d) + @b)                      # λ = W+ (A^T d + b)
      @xs      = @h_inv.dot(@a.dot(lambdas) - @d)                       # x = H^(-1) (A λ - d)

      return @xs, @lambdas
    end

    private
    def constraint_matrix
      a = []
      b = []
      @constraint.each do |c|
        a << c.a.to_flat_array
        b << c.b
      end
      @a = NMatrix.new [a.size, @size], a.flatten
      @b = NMatrix.new [@size, 1]
      return @a, @b
    end
  end

  ##
  # Class that represents a general quadratic optimizer with active set search.
  class GeneralQuadraticOptimizer < Algorithm
    # Active set
    attr_reader :active_set

    ##
    # Initializer for our quadratic optimizer with active set
    def initialize(options, objective)
      raise ArgumentError, "Objective function must be a QuadraticObjectiveFunction" unless objective.is_a? QuadraticObjectiveFunction
      super options, objective
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
      raise ArgumentError, "Initial guess must be a vector" unless x0.is_a? NMatrix
      raise ArgumentError, "x0 has a wrong size. Must be of size #{@size} x 1" unless (x0.shape == [@size, 1])
      binding.pry if ARGV.include? "-d"
      puts "Starting optimization problem".green if @options[:debug] > 0
      reorder_constraints

      puts " - Cholesky factorization for H matrix and g vector".green if @options[:debug] > 0
      preparation

      x = x0

      puts " - Evaluating initial active set".green if @options[:debug] > 0
      active_set(x)

      puts " - Starting solution loop".green if @options[:debug] > 0
      iter = 0

      while (iter != @options[:iterations]) do
        iter += 1

        puts " -- Soving quadratic programming problem".green if @options[:debug] > 1
        z, lambdas = q_problem
        puts "    z = #{x.to_flat_array}".yellow if @options[:debug] > 1
        puts "    λ = #{(lambdas ? lambdas.to_flat_array : 'empty vector')}".yellow if @options[:debug] > 1

        if not feasible?(z)
          puts "  --- Point is not feasible.".green if @options[:debug] > 1
          x = line_problem(z, x)
          #active_set(x)

        else
          x = z
          feasible_set(x, lambdas)
          if @feasible_set == []
            break
          else
            #active_set(x)
            @active_set -= @feasible_set
          end
        end
      end
      binding.pry if ARGV.include? "-d"
      return x
    end

    private
    ##
    # Returns solution for the sub-problem
    def q_problem
      if @active_set != []
        a, b = active_set_matrix
        at   = a.transpose.dup

        k    = at.dot (@hinv.dot a)
        kinv = k.invert

        lambdas = kinv.dot ((at.dot @d) - b)
        xs      = (@hinv.dot (a.dot(lambdas))) - @d
      else
        lambdas = nil
        xs = -@d
      end
      binding.pry if ARGV.include? "-d"
      return xs, lambdas
    end

    ##
    # Executes a line search to find the best point on the contraints
    def line_problem(z, x)
      puts " ---- Performing line search".green if @options[:debug] > 2
      puts "      x = #{x.to_flat_array}".yellow if @options[:debug] > 2
      puts "      z = #{z.to_flat_array}".yellow if @options[:debug] > 2
      t_set = []
      @constraints.each do |c|
        puts " ---- New constraint checked" if @options[:debug] > 4
        if (c.type == :inequality and c.f(z) < 0.0)
          print " ---- Constraint violated: ".green if @options[:debug] > 2
          t = (-c.b - (c.a.transpose.dot(x))[0])/((c.a.transpose.dot((x - z)))[0])
          #raise RuntimeError, "t is too big: 0 <= t = #{t} <= 1" if (t < 0 or t > 1)
          t_set << t
          puts " t = #{t}".yellow if @options[:debug] > 2
        end
      end
      t_set.each_with_index { |t,i|; (@active_set - [i]) + [i] if t == t_set.min }
      ret = x + ((x - z) * t_set.min)
      puts " ---- New point considered: ".green if @options[:debug] > 2
      puts "      x = #{ret.to_flat_array}".yellow if @options[:debug] > 2
      return ret
    end

    ##
    # Evaluates the active set for the defined step. Also evaluate active inequality set
    def active_set(x)
      puts " ---- Evaluating active set".green if @options[:debug] > 2
      @active_set = []
      @constraints.each_with_index do |c, j|
        case c.type
          when :equality
            @active_set << j
          when :inequality
            @active_set << j if c.active? x
        end
      end
      puts "     A = #{@active_set}".yellow if @options[:debug] > 0
    end

    ##
    # Creates a new feasible set
    def feasible_set(x, lambdas)
      puts " ---- Checking feasible set".green if @options[:debug] > 4
      @feasible_set = []
      @active_set.each_with_index do |j, i|
        @feasible_set << j if lambdas[i] < 0 # -@options[:tolerance]
      end
      puts "      Φ = #{@feasible_set}".yellow if @options[:debug] > 4
    end

    ##
    # Extracts matrix A and vector b given the active set
    def active_set_matrix
      a_ary = []
      b_ary = []
      @active_set.each do |j|
        a_ary << @constraints[j].a.to_flat_array
        b_ary << @constraints[j].b
      end
      a_matrix = NMatrix.new([@size, @active_set.size], a_ary.flatten)
      b_matrix = NMatrix.new([@active_set.size, 1], b_ary.flatten)
      return a_matrix, b_matrix
    end

    ##
    # Evaluate Cholesky decomposition, used in subproblem solution
    def preparation
      @hinv = @objective.s.invert
      @g = @objective.b
      @d = @hinv.dot(@g)
    end

    ##
    # Move to the first positions equality constraints. Used to define the active set.
    def reorder_constraints
      puts " - Reordering constraints".green if @options[:debug] > 0
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
  end
end
