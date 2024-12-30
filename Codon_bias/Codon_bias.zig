const std = @import("std");
const print = std.debug.print;

const Sequence = struct{
    name : []const u8,
    letters : []const u8,
    current: usize,
    pub fn next(Self: *Sequence) ?u8{
        if(Self.current >= Self.letters.len) return null;
        if(Self.letters[Self.current] == '\n') Self.current += 1; //skip newlines
        if(Self.letters[Self.current] == '\r') Self.current += 1; //skip carriage returns
        if(Self.current >= Self.letters.len) return null; //check again
        defer Self.current += 1;
        return switch(Self.letters[Self.current]){
            'U' => 'T', //return timin instead of uracil
            'u' => 't',
            'N' => '-', //return - instead of other "missing" markers
            'n' => '-',
            'X' => '-',
            'x' => '-',
            else=>Self.letters[Self.current],
        };
    }
    pub fn jump(Self: *Sequence, n: usize) void{
        var count: usize = 0;
        while(count < n){
            count += 1;
            Self.current += 1;
            if(Self.letters[Self.current] == '\n') Self.current += 1;
            if(Self.letters[Self.current] == '\r') Self.current += 1;
        }
    }
    pub fn nextBatch(Self: *Sequence, comptime n:usize, target: *[n]u8) void{
        for(target, 0..) |t, i|{
            _ = t;
            target[i] = Self.next() orelse '-';
        }
    }
    pub fn nextCodon(Self: *Sequence) ?[3]u8{
        var resp = [_]u8{0,0,0};
        var letter: ?u8 = Self.next();
        if(letter == null) return null;
        resp[0] = letter.?;
        letter = Self.next();
        if(letter == null) return null;
        resp[1] = letter.?;
        letter = Self.next();
        if(letter == null) return null;
        resp[2] = letter.?;
        return resp;
    }
    pub fn start(Self: *Sequence) void{
        Self.current = 0;
    }
    //only for testing because my test data has stupid long names
    pub fn shortName(Self: *const Sequence) []const u8{
        for(Self.name,0..) |l, i|{
            if(l == '[') return Self.name[0..i];
        }
        return Self.name;
    }
};

const CodonTable = struct{
    codons: []u8,
    amino_acids: []u8,
    ///returns the amino acid corresponding to the given codon(insert bit magic here)
    fn translate(self: *CodonTable, codon: u8) u8{
        for(self.codons, self.amino_acids) |c, a|{
            if(c == codon) return a;
        }
        return 'X';
    }
    ///returns a slice from the codons containing all possible encodings for the given amino acid
    fn getCodons(self: *CodonTable, amino_acid: u8) []u8{
        var start:usize = 0;
        while(self.amino_acids[start] != amino_acid){start += 1;} //find beginning of codons
        var end:usize = start+1;
        while(self.amino_acids[end] != amino_acid){end += 1;} //find end of codons
        return self.codons[start..end];
    }
};
pub fn main() !void{
    //allocator
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    
    //arg reading
    var args = try std.process.argsWithAllocator(allocator);
    _=args.skip();
    var fname_nucleotide: ?[]const u8 = null;
    var tname: ?[]const u8 = null;
    var fname_output: ?[]const u8 = null;
    while(args.next()) |arg|{
        if(std.mem.eql(u8, arg, "-n")){
            fname_nucleotide = args.next();
        }
        else if(std.mem.eql(u8, arg, "-t")){
            tname = args.next();
        }
        else if(std.mem.eql(u8, arg, "-h")){
            print(
            \\switch syntax: -<switch> <parameter>
            \\-h
            \\  displays this message
            \\-n <file name>
            \\  specifies the file containing the sequences as nucleotides
            \\-t <table name>
            \\  choose witch translation table to use for codons
            \\-o <file name>
            \\  modifies the name of the file(with extension .csv) to write the results in if the file already exists its contents will be overwritten
            \\  the default name is: codon_bias.csv
            \\  write "none", if you don't want to write the contents to a file
            \\
            , .{});
            return;
        }
        else if(std.mem.eql(u8, arg, "-o")){
            fname_output = args.next();
        }
        //std.debug.print("{s}\n",.{arg});
    }
    if(fname_nucleotide == null){
        print("bruh", .{});
        return;
    }
}

//note: no toUppercase, just if(c > 'Z') c -= 32;
fn packCodons(source: *Sequence, target: *Sequence) void{
    var i: usize = 0;
    while(source.nextCodon()) |codon|{
        if(codon[0] == '-') target[i] = 'X';
        
        i += 1;
    }
}
fn fillSequence(target: []Sequence, source: []const u8) usize{
    var count: usize = 0;
    var start: usize = 0;
    var i: usize = 0;
    
    while(i < source.len) : (count += 1){
        if(count == 100){
            std.debug.print("!!---only first 100 sequences read---!!", .{});
            break;
        }
        //name
        i += 1; //skip >
        start = i;
        while(source[i] != '\n'){i+=1;}//skip to new line
        target[count].name = source[start..i];
        i += 1; //go over new line
        if(source[i] == '\r') i+=1; //ignore carriage return (fuck u microsoft)
        //sequence
        start = i;
        while(i < source.len){
            if(source[i] == '>') break;
            i+=1;
        }
        //handle extra characters
        const extra_spaces: u8 = switch(source[i-1]){
            '\r' => 2,
            '\n' => 1,
            else => 0,
        };
        target[count].letters = source[start..(i-extra_spaces)];
    }
    return count;
}
fn readFile(filename: []const u8, allocator: std.mem.Allocator) ![]u8{
    const cwd = std.fs.cwd();
    var file = try cwd.openFile(filename, .{});
    defer file.close();
    const stat = try file.stat();
    const size = stat.size;
    const content = try allocator.alloc(u8, size);
    _ = try file.read(content);
    return content;
}